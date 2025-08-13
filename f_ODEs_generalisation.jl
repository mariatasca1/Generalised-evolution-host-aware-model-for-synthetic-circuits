function f_ODEs_generalisation(du, y, parameters_matlab, t)
   
    params = parameters_matlab[1];
    M = parameters_matlab[2];
    subpops_up = parameters_matlab[3];
    z_values_up = parameters_matlab[4];
    z_values_dn = parameters_matlab[5];
    # Unpack parameters

    (n, num_genes, N, mass, nut, nq, v_e, K_e, n_R, n_Z, n_Q) = params[1:11]; 
    n = Int(n);
    num_genes = Int(num_genes);
    m = 1 + 8*n # to determine the size of G
    G = similar(y, m, num_genes); # G is the matrix that will have all the parameters needed for genes

    G[1,:] = collect(params[12:12+num_genes-1]); # this is to get the n_hA... on the first row

    (v_TX_R, v_TX_Z, v_TX_Q, K_TX_R, K_TX_nR, K_Q, h_Q, m_deg, kb_TL, ku_TL, v_TL, K_TL, K_H, h_H, p_deg) = params[12+num_genes:27+num_genes-1]; #next set of parameters

    x = 27 + num_genes; #to make indices tracking easier
    matrix_param = params[x];
    #@show typeof(matrix_param)
    #@show size(matrix_param)


    for i in 1:num_genes
        G[ 2:n+1, i] = matrix_param[2*i-1, :]; #  will help with transcription rates later, alpha coefficients

        G[ n+2:2*n+1, i] = matrix_param[2*i, :]; # beta coefficients, will help in the m_X equations
        
    end
    

   # Now onto variables!

   x = 1  # Initialize for automatic numbering
   subpop = y[x:x+n-1]; x += n;
   e = y[x:x+n-1]; x += n;
   m_R = y[x:x+n-1]; x += n;
   m_Z = y[x:x+n-1]; x += n;
   m_Q = y[x:x+n-1]; x += n;

   for i in 1: num_genes
       G[2*n+2:3*n+1, i] = y[x:x+n-1]; #we introduce the m_HA, m_HB, m_HC... as column vectors into G, but when used we will transpose them back
       x = x + n;
   end

   TL_R = y[x:x+n-1]; x=x+n;
   TL_Z = y[x:x+n-1]; x=x+n;
   TL_Q = y[x:x+n-1]; x=x+n;

   for i in 1:num_genes
    G[3*n+2:4*n+1, i] = y[x:x+n-1]; #we introduce the TL_HA, TH_HB, TH_HC... as column vectors into G, but when used we will transpose them back
    x = x + n;
   end

   R = y[x:x+n-1]; x=x+n;
   Z = y[x:x+n-1]; x=x+n;
   Q = y[x:x+n-1]; x=x+n;

   for i in 1:num_genes
    G[4*n+2:5*n+1, i] = y[x:x+n-1]; #we introduce the HA, HB, HC ... as column vectors into G, but when used we will transpose them back
    x = x + n;
   end

    #Calculate rates

    #Energy
    epsilon = nq*Z*v_e*nut ./ (K_e + nut);
    #@show epsilon

    #Autoregulation of m_Q
    IQ = (1 ./(1 .+ (Q/K_Q).^h_Q));   #Hill function
    #@show IQ
    
    # Transcription regulation for synthetic genes. Choose one set to 
    # uncomment, or make a new set that suits whatever regulatory network
    # topology you want to implement.

    # IHA = 1 no repression (ig a vector of 1s?)

    # M is a matrix that records all interactions between the components of a circuit
    # i.e. if gene 2 inhibits gene 3 then element [3,2] = -1
    # if gene 2 promotes gene 3 then element [3,2] = 1
    # a repressilator would have the matrix [(0, 0, -1), (-1, 0, 0), (0, -1, 0)]

    for i in 1:num_genes
            for j in 1:num_genes
                if M[i,j] == -1
                    G[5*n+2:6*n+1, i] = (1 ./ (1 .+ (G[4*n+2:5*n+1, j] /K_H).^h_H));
                else
                    if M[i,j] == 1
                        G[5*n+2:6*n+1, i] = 1 .-(1 ./ (1 .+ (G[4*n+2:5*n+1, j] /K_H).^h_H));
                    end
                end
            end
    end

    # Transcription rates

    w_R = (v_TX_R * e) ./ (K_TX_R .+ e);
    #@show w_R
    w_Z = (v_TX_Z * e) ./ (K_TX_nR .+ e);
    #@show w_Z
    w_Q = (v_TX_Q * e) ./ (K_TX_nR .+ e) .* IQ;
    #@show w_Q

    for i in 1:num_genes
        G[6*n+2:7*n+1, i] = (G[2:n+1,i] .*e) ./ (K_TX_nR .+ e) .* G[5*n+2:6*n+1, i];
    end

    # Translation

    TL_rate = (v_TL * e) ./ (K_TL .+ e);
    #@show TL_rate
    
    gamma_R = TL_R .* TL_rate / n_R;
    #@show gamma_R

    gamma_Z = TL_Z .* TL_rate / n_Z;
    #@show gamma_Z

    gamma_Q = TL_Q .* TL_rate / n_Q;
    #@show gamma_Q

    for i in 1:num_genes
        G[7*n+2:8*n+1, i] = G[3*n+2: 4*n + 1, i] .* TL_rate ./ G[1, i];
    end
    # @show G
  

    # Useful summations
    m_all = m_R + m_Z + m_Q;
    for i in 1:num_genes
        m_all = m_all + G[2*n+2:3*n+1, i];
    end
    # @show m_all

    TL_all = TL_R + TL_Z + TL_Q;
    for i in 1:num_genes
        TL_all = TL_all .+ G[3*n+2:4*n+1, i];
    end
    # @show TL_all

    gamma_all = gamma_R + gamma_Z + gamma_Q;
    for i in 1:num_genes
        gamma_all = gamma_all .+ G[7*n+2:8*n+1, i];
    end
    # show gamma_all

    # Growth rate
    GR = TL_rate .* TL_all / mass;
    # @show GR
    # Dilution. To ensure that the numerical implementation follows the
    # mathematical formulation from Equation (1), it can be written without
    # the if/else structure. If using the if/else structure (i.e. the
    # mathematically correct formulation), then abs/rel tolerances usually
    # need to be lowered for the simulations to work.

    buffer = sum(subpop) - N;
    # @show buffer
    # if sum(subpop) > N
    #     buffer = sum(subpop) - N;
    # else
    #     buffer = 0;
    # end

    # ODEs

    # For explanations of the 'number of cells in each state', see
    # Section S1.2. For the cell variables, see Section S1.1.

    dydt = similar(y); # Initialise to store values
    x = 1; # Initialise for automatic numbering
    subpops_up_int = [Int.(v) for v in subpops_up]

    #@show subpops_up
    #@show typeof(subpops_up)


    # Number of cells in each state
    for i in 1:n
        idxs = subpops_up_int[i]

             # First, check if idxs is not empty to avoid errors
        if !isempty(idxs)
               # Explicit elementwise multiplication with a comprehension, multiplying each element
               
            term1 = sum(subpop[idxs] .* GR[idxs] .* z_values_up[i])
            
        else
            term1 = zero(eltype(subpop))  # zero dual number of the right type
         end

        term2 = subpop[i] * GR[i] * (1 - sum(z_values_dn[i]))  # Normal self-division
        term3 = subpop[i] * buffer  # Dilution
       

        dydt[i] = term1 + term2 - term3
    end
    x=x+n;
    #@show dydt[1:8]
    

    dydt[x:x+n-1] = epsilon - (gamma_R*n_R + gamma_Z*n_Z + gamma_Q*n_Q) - GR.*e; #@show dydt[x:x+n-1]
    for i in 1:num_genes
        dydt[x:x+n-1] = dydt[x:x+n-1] - G[7*n+2:8*n+1, i] * G[1,i];
        #@show dydt[x:x+n-1]
        #@show G[7*n+2:8*n+1,i] * G[1,i]
    end # the rest of e
    #@show dydt[x:x+n-1]
    x = x+n;
    #@show epsilon
    #@show gamma_R * n_R
    #@show gamma_Z *n_Z
    #@show gamma_Q *n_Q
    #@show GR .* e
 

    # m_R
    dydt[x:x+n-1] = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR.+m_deg).*m_R; 
    x=x+n;

    # m_Z
    dydt[x:x+n-1] = w_Z - kb_TL*m_Z.*R + ku_TL*TL_Z + gamma_Z - (GR.+m_deg).*m_Z; 
    x=x+n;

    # m_Q
    dydt[x:x+n-1] = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR.+m_deg).*m_Q; 
    x=x+n;

    #m_HAs:
    for i in 1:num_genes
        dydt[x:x+n-1] = G[6*n+2: 7*n+1, i] - G[n+2:2*n+1, i] .* G[2*n+2: 3*n+1, i] .* R + ku_TL * G[3*n+2: 4*n+1, i] + G[7*n+2:8*n+1, i] - (GR .+ m_deg) .* G[2*n+2:3*n+1, i];
        x = x+n;
    end

    # TL_R
    dydt[x:x+n-1] = kb_TL*m_R.*R - ku_TL*TL_R - gamma_R - GR.*TL_R; x=x+n;
    
    # TL_Z
    dydt[x:x+n-1] = kb_TL*m_Z.*R - ku_TL*TL_Z - gamma_Z - GR.*TL_Z; x=x+n;
    
    # TL_Q
    dydt[x:x+n-1] = kb_TL*m_Q.*R - ku_TL*TL_Q - gamma_Q - GR.*TL_Q; x=x+n;
    
    #TL_HAs:
    for i in 1:num_genes
        dydt[x:x+n-1] = G[n+2:2*n+1, i] .*G[2*n+2:3*n+1, i] .*R - ku_TL* G[3*n+2:4*n+1, i] - G[7*n+2:8*n+1, i] - GR.* G[3*n+2:4*n+1, i];
        x = x+n;
    end
    
    # R
    #@show x
    #@show gamma_R, ku_TL, TL_all, kb_TL, m_all, R, gamma_all, GR
    dydt[x:x+n-1] = gamma_R + ku_TL* (TL_Q +TL_R + TL_Z) - kb_TL*(m_R + m_Q + m_Z).*R + gamma_all - GR.*R;
    #@show dydt[x:x+n-1]
    #@show G[n+2:2*n+1, :],G[2*n+2:3*n+1, :], G[3*n+2:4*n+1, :]

    for i in 1:num_genes
        dydt[x:x+n-1] = dydt[x:x+n-1] + ku_TL * G[3*n+2:4*n+1, i] - G[n+2:2*n+1, i].* G[2*n+2: 3*n+1,i].*R
        #@show i
        #@show ku_TL * G[3*n+2:4*n+1, i] - G[n+2:2*n+1, i].* G[2*n+2: 3*n+1].*R
        #@show ku_TL * G[3*n+2:4*n+1, i]
        #@show G[n+2:2*n+1, i].* G[2*n+2: 3*n+1, i].*R
        #@show G[n+2:2*n+1, i]
        #@show G[2*n+2:3*n+1,i]
        #@show R
        
    end #the rest of R
    x = x+n;

    # Z
    dydt[x:x+n-1] = gamma_Z - GR.*Z; x=x+n;
    #@show dydt[x:x+n-1]
    #@show typeof(dydt[x:x+n-1])

    # Q
    # @show gamma_Q, Q
    dydt[x:x+n-1] = gamma_Q - GR.*Q; 
    x=x+n;
    #@show dydt[x:x+n-1]
    #@show typeof(dydt[x:x+n-1])

    #HAs: 
    #@show GR, p_deg
    #@show G[7*n+2:8*n+1,:]
    #@show G[4*n+2:5*n+1,:]

    for i in 1:num_genes
        dydt[x:x+n-1] = G[7*n+2:8*n+1, i] - (GR .+ p_deg) .* G[4*n+2:5*n+1, i]; 
        x = x+n;
    end

    #@show dydt[x:x+n-1]
    #@show typeof(dydt[x:x+n-1])


    #@show length(dydt)
    #@show typeof(dydt)
    #@show dydt

    du[:] = dydt
    # @show length(du)
    #@show typeof(du)
    #@show du
    
end