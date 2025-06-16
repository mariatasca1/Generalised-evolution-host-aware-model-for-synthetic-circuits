function f_ODEs_generalisation(t,y, params, M, subpops_up,z_values_up, z_values_dn)
    @info "t = $t"
    # Unpack parameters
    (n, num_genes, N, mass, nut, nq, v_e, K_e, n_R, n_Z, n_Q) = params[1:11]; 
    m = 1 + 8*n # to determine the size of G
    G = zeros(m, Int(num_genes)); # G is the matrix that will have all the parameters needed for genes

    G[1,:] = collect(params[12:12+num_genes-1]); # this is to get the n_hA... on the first row

    (v_TX_R, v_TX_Z, v_TX_Q, K_TX_R, K_TX_nR, K_Q, h_Q, m_deg, kb_TL, ku_TL, v_TL, K_TL, K_H, h_H, p_deg) = params[12+num_genes:27+num_genes-1]; #next set of parameters

    x = 27 + num_genes; #to make indices tracking easier

    for i in 1:Int(num_genes)
        G[ 2:n+1, i] = params[x][2*i-1, :]; # Must be column vectors, will help with transcription rates later, alpha coefficients

        G[ n+2:2*n+1, i] = params[x][2*i, :]; # beta coefficients, will help in the m_X equations
        
    end

    # Now onto variables!

    x = 1  # Initialize for automatic numbering
   subpop = y[x:x+n-1]; x += n;
   e = y[x:x+n-1]; x += n;
   m_R = y[x:x+n-1]; x += n;
   m_Z = y[x:x+n-1]; x += n;
   m_Q = y[x:x+n-1]; x += n;

   for i in 1: Int(num_genes)
       G[2*n+2:3*n+1, i] = y[x:x+n-1]; #we introduce the m_HA, m_HB, m_HC... as column vectors into G, but when used we will transpose them back
       x = x + n;
   end

   TL_R = y[x:x+n-1]; x=x+n;
   TL_Z = y[x:x+n-1]; x=x+n;
   TL_Q = y[x:x+n-1]; x=x+n;

   for i in 1:Int(num_genes)
    G[3*n+2:4*n+1, i] = y[x:x+n-1]; #we introduce the TL_HA, TH_HB, TH_HC... as column vectors into G, but when used we will transpose them back
    x = x + n;
   end

   R = y[x:x+n-1]; x=x+n;
   Z = y[x:x+n-1]; x=x+n;
   Q = y[x:x+n-1]; x=x+n;

   for i in 1:Int(num_genes)
    G[4*n+2:5*n+1, i] = y[x:x+n-1]; #we introduce the HA, HB, HC ... as column vectors into G, but when used we will transpose them back
    x = x + n;
   end

    #Calculate rates

    #Energy
    epsilon = nq*Z*v_e*nut ./ (K_e + nut);

    #Autoregulation of m_Q
    IQ = (1 ./(1 .+ (Q/K_Q).^h_Q));   #Hill function
  
    # Transcription regulation for synthetic genes. Choose one set to 
    # uncomment, or make a new set that suits whatever regulatory network
    # topology you want to implement.

    # IHA = 1 no repression (ig a vector of 1s?)

    # M is a matrix that records all interactions between the components of a circuit
    # i.e. if gene 2 inhibits gene 3 then element [3,2] = -1
    # if gene 2 promotes gene 3 then element [3,2] = 1
    # a repressilator would have the matrix [(0, 0, -1), (-1, 0, 0), (0, -1, 0)]
    for i in 1:Int(num_genes)
            for j in 1:Int(num_genes)
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
    w_Z = (v_TX_Z * e) ./ (K_TX_nR .+ e);
    w_Q = (v_TX_Q * e) ./ (K_TX_nR .+ e) .* IQ;

    for i in 1:Int(num_genes)
        G[6*n+2:7*n+1, i] = (G[2:n+1,i] .*e) ./ (K_TX_nR .+ e) .* G[5*n+2:6*n+1, i];
    end

    # Translation
    TL_rate = (v_TL * e) ./ (K_TL .+ e);
    @info "e = $e"
    gamma_R = TL_R .* TL_rate / n_R;
    gamma_Z = TL_Z .* TL_rate / n_Z;
    gamma_Q = TL_Q .* TL_rate / n_Q;

    for i in 1:Int(num_genes)
        G[7*n+2:8*n+1, i] = G[3*n+2: 4*n + 1, i] .* TL_rate ./ G[1, i];
    end

    # Useful summations
    m_all = m_R + m_Z + m_Q;
    for i in 1:Int(num_genes)
        m_all = m_all + G[2*n+2:3*n+1, i];
    end

    TL_all = TL_R + TL_Z + TL_Q;
    for i in 1:Int(num_genes)
        TL_all = TL_all .+ G[3*n+2:4*n+1, i];
    end

    gamma_all = gamma_R + gamma_Z + gamma_Q;
    for i in 1:Int(num_genes)
        gamma_all = gamma_all .+ G[7*n+2:8*n+1, i];
    end

    # Growth rate
    GR = TL_rate .* TL_all / mass;
    # Dilution. To ensure that the numerical implementation follows the
    # mathematical formulation from Equation (1), it can be written without
    # the if/else structure. If using the if/else structure (i.e. the
    # mathematically correct formulation), then abs/rel tolerances usually
    # need to be lowered for the simulations to work.

    buffer = sum(subpop) - N;
    # if sum(subpop) > N
    #     buffer = sum(subpop) - N;
    # else
    #     buffer = 0;
    # end

    # ODEs

    # For explanations of the 'number of cells in each state', see
    # Section S1.2. For the cell variables, see Section S1.1.

    dydt = zeros(Int(length(y))); # Initialise to store values
    x = 1; # Initialise for automatic numbering

    # Number of cells in each state
    for i in 1:n
        idxs = Int.(subpops_up[i]);
        dydt[i] = sum(subpop[idxs].*GR[idxs].*z_values_up[i])  # Division from upstream states
            + subpop[i]*GR[i]*(1-sum(z_values_dn[i]))  # Normal 'self-division'
            - subpop[i]*buffer; # Dilution
    end
    x=x+n;
    
    # e
    dydt[x:x+n-1] = epsilon - (gamma_R*n_R + gamma_Z*n_Z + gamma_Q*n_Q)
              - GR.*e;
    for i in 1:Int(num_genes)
        dydt[x:x+n-1] = dydt[x:x+n-1] + G[7*n+2:8*n+1, i] .* G[1,i];
    end
    x = x+n;
    # m_R
    dydt[x:x+n-1] = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR.+m_deg).*m_R; x=x+n;
    # m_Z
    dydt[x:x+n-1] = w_Z - kb_TL*m_Z.*R + ku_TL*TL_Z + gamma_Z - (GR.+m_deg).*m_Z; x=x+n;
    # m_Q
    dydt[x:x+n-1] = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR.+m_deg).*m_Q; x=x+n;

    for i in 1:Int(num_genes)
        dydt[x:x+n-1] = G[6*n+2: 7*n+1, i] - G[n+2:2*n+1, i] .* G[2*n+2: 3*n+1, i] .* R + ku_TL * G[3*n+2: 4*n+1, i] 
        + G[7*n+2:8*n+1, i] - (GR .+ m_deg) .* G[2*n+2:3*n+1, i];
        x = x+n;
    end
    
    # TL_R
    dydt[x:x+n-1] = kb_TL*m_R.*R - ku_TL*TL_R - gamma_R - GR.*TL_R; x=x+n;
    # TL_Z
    dydt[x:x+n-1] = kb_TL*m_Z.*R - ku_TL*TL_Z - gamma_Z - GR.*TL_Z; x=x+n;
    # TL_Q
    dydt[x:x+n-1] = kb_TL*m_Q.*R - ku_TL*TL_Q - gamma_Q - GR.*TL_Q; x=x+n;

    for i in 1:Int(num_genes)
        dydt[x:x+n-1] = G[n+2:2*n+1, i] .*G[2*n+2:3*n+1, i] .*R - ku_TL* G[3*n+2:4*n+1, i] - G[7*n+2:8*n+1, i] 
        - GR.* G[3*n+2:4*n+1, i];
        x = x+n;
    end

    # R
    dydt[x:x+n-1] = gamma_R + ku_TL*TL_all - kb_TL*m_all.*R + gamma_all - GR.*R;

    for i in 1:Int(num_genes)
        dydt[x:x+n-1] = dydt[x:x+n-1] - ku_TL * G[3*n+2:4*n+1, i] + kb_TL * G[2*n+2:3*n+1, i] .*R 
                      + ku_TL * G[3*n+2:4*n+1, i] - G[n+2:2*n+1, i].* G[2*n+2: 3*n+1].*R
    end
    x = x+n;

    # Z
    dydt[x:x+n-1] = gamma_Z - GR.*Z; x=x+n;
    # Q
    dydt[x:x+n-1] = gamma_Q - GR.*Q; x=x+n;

    for i in 1:Int(num_genes)
        dydt[x:x+n-1] = G[7*n+2:8*n+1, i] - (GR .+ p_deg) .* G[5*n+2:6*n+1, i] 
    end

    return dydt
end