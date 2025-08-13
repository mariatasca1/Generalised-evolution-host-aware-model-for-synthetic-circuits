function f_ODEs_3genes_original(du, y, parameters_matlab, t)
    # Unpack parameters
    params = parameters_matlab[1];
    subpops_up = parameters_matlab[2];
    z_values_up = parameters_matlab[3];
    z_values_dn = parameters_matlab[4];
    x = 1; # Initialise for automatic numbering
    n = Int(params[x]); x=x+1;
    N = params[x]; x=x+1;
    mass = params[x]; x=x+1;
    nut = params[x]; x=x+1;
    nq = params[x]; x=x+1;
    v_e = params[x]; x=x+1;
    K_e = params[x]; x=x+1;
    n_R = params[x]; x=x+1;
    n_Z = params[x]; x=x+1;
    n_Q = params[x]; x=x+1;
    n_HA = params[x]; x=x+1;
    n_HB = params[x]; x=x+1;
    n_HC = params[x]; x=x+1;
    v_TX_R = params[x]; x=x+1;
    v_TX_Z = params[x]; x=x+1;
    v_TX_Q = params[x]; x=x+1;
    K_TX_R = params[x]; x=x+1;
    K_TX_nR = params[x]; x=x+1;
    K_Q = params[x]; x=x+1;
    h_Q = params[x]; x=x+1;
    m_deg = params[x]; x=x+1;
    kb_TL = params[x]; x=x+1;
    ku_TL = params[x]; x=x+1;
    v_TL = params[x]; x=x+1;
    K_TL = params[x]; x=x+1;
    K_H = params[x]; x=x+1;
    h_H = params[x]; x=x+1;
    p_deg = params[x]; x=x+1;
    alpha_A_vec = params[x:x+n-1]; x=x+length(alpha_A_vec); 
    beta_A_vec = params[x:x+n-1]; x=x+length(beta_A_vec);
    alpha_B_vec = params[x:x+n-1]; x=x+length(alpha_B_vec);
    beta_B_vec = params[x:x+n-1]; x=x+length(beta_B_vec);
    alpha_C_vec = params[x:x+n-1]; x=x+length(alpha_C_vec);
    beta_C_vec = params[x:x+n-1];

    # Now onto variables!

   x = 1  # Initialize for automatic numbering
   subpop = y[x:x+n-1]; x += n;
   e = y[x:x+n-1]; x += n;
   m_R = y[x:x+n-1]; x += n;
   m_Z = y[x:x+n-1]; x += n;
   m_Q = y[x:x+n-1]; x += n;

    m_HA = y[x:x+n-1]; x=x+n;
    m_HB = y[x:x+n-1]; x=x+n;
    m_HC = y[x:x+n-1]; x=x+n;

   TL_R = y[x:x+n-1]; x=x+n;
   TL_Z = y[x:x+n-1]; x=x+n;
   TL_Q = y[x:x+n-1]; x=x+n;

    TL_HA = y[x:x+n-1]; x=x+n;
    TL_HB = y[x:x+n-1]; x=x+n;
    TL_HC = y[x:x+n-1]; x=x+n;

   R = y[x:x+n-1]; x=x+n;
   Z = y[x:x+n-1]; x=x+n;
   Q = y[x:x+n-1]; x=x+n;

    HA = y[x:x+n-1]; x=x+n;
    HB = y[x:x+n-1]; x=x+n;
    HC = y[x:x+n-1];

    #Calculate rates

    #Energy
    epsilon = nq*Z*v_e*nut ./ (K_e + nut);
   
    #Autoregulation of m_Q
    IQ = (1 ./(1 .+ (Q/K_Q).^h_Q));   #Hill function
  
    # Transcription regulation for synthetic genes. Choose one set to 
    # uncomment, or make a new set that suits whatever regulatory network
    # topology you want to implement.

    # IHA = 1 no repression (ig a vector of 1s?)

  
    IHA = (1 ./ (1 .+ (HC/K_H).^h_H)); # Repression of m_HA from HC
    IHB = (1 ./ (1 .+ (HA/K_H).^h_H)); # Repression of m_HB from HA
    IHC = (1 ./ (1 .+ (HB/K_H).^h_H)); # Repression of m_HC from HB
   

    # Transcription rates
    w_R = (v_TX_R * e) ./ (K_TX_R .+ e);
    w_Z = (v_TX_Z * e) ./ (K_TX_nR .+ e);
    w_Q = (v_TX_Q * e) ./ (K_TX_nR .+ e) .* IQ;
 

    w_HA = (alpha_A_vec .* e) ./ (K_TX_nR .+ e) .* IHA;
    w_HB = (alpha_B_vec .* e) ./ (K_TX_nR .+ e) .* IHB;
    w_HC = (alpha_C_vec .* e) ./ (K_TX_nR .+ e) .* IHC;
   
    # Translation
    TL_rate = (v_TL * e) ./ (K_TL .+ e);
    gamma_R = TL_R .* TL_rate / n_R;
    gamma_Z = TL_Z .* TL_rate / n_Z;
    gamma_Q = TL_Q .* TL_rate / n_Q;
  

    gamma_HA = TL_HA .* TL_rate ./ n_HA;
    gamma_HB = TL_HB .* TL_rate ./ n_HB;
    gamma_HC = TL_HC .* TL_rate ./ n_HC;
   

    # Useful summations
    m_all = m_R + m_Z + m_Q + m_HA + m_HB + m_HC;

    TL_all = TL_R + TL_Z + TL_Q + TL_HA + TL_HB + TL_HC;

    gamma_all = gamma_R + gamma_Z + gamma_Q + gamma_HA + gamma_HB + gamma_HC;

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

    dydt = similar(y); # Initialise to store values
    x = 1; # Initialise for automatic numbering
    subpops_up_int = [Int.(v) for v in subpops_up]

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
    # e
    dydt[x:x+n-1] = epsilon - (gamma_R*n_R + gamma_Z*n_Z + gamma_Q*n_Q + gamma_HA .*n_HA + gamma_HB .* n_HB + gamma_HC .* n_HC)
              - GR.*e;x = x + n;
   
    # m_R
    dydt[x:x+n-1] = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR.+m_deg).*m_R; x=x+n;
    # m_Z
    dydt[x:x+n-1] = w_Z - kb_TL*m_Z.*R + ku_TL*TL_Z + gamma_Z - (GR.+m_deg).*m_Z; x=x+n;
    # m_Q
    dydt[x:x+n-1] = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR.+m_deg).*m_Q; x=x+n;

    # m_HA
    dydt[x:x+n-1] = w_HA .- beta_A_vec.*m_HA.*R .+ ku_TL*TL_HA .+ gamma_HA .- (GR.+m_deg).*m_HA; x=x+n;
    #  m_HB
    dydt[x:x+n-1] = w_HB .- beta_B_vec.*m_HB.*R .+ ku_TL*TL_HB .+ gamma_HB .- (GR.+m_deg).*m_HB; x=x+n;
    # m_HC
    dydt[x:x+n-1] = w_HC .- beta_C_vec.*m_HC.*R .+ ku_TL*TL_HC .+ gamma_HC .- (GR.+m_deg).*m_HC; x=x+n;
    
    # TL_R
    dydt[x:x+n-1] = kb_TL*m_R.*R .- ku_TL*TL_R .- gamma_R .- GR.*TL_R; x=x+n;
    # TL_Z
    dydt[x:x+n-1] = kb_TL*m_Z.*R .- ku_TL*TL_Z .- gamma_Z .- GR.*TL_Z; x=x+n;
    # TL_Q
    dydt[x:x+n-1] = kb_TL*m_Q.*R .- ku_TL*TL_Q .- gamma_Q .- GR.*TL_Q; x=x+n;
    
    # TL_HA
    dydt[x:x+n-1] = beta_A_vec.*m_HA.*R .- ku_TL*TL_HA .- gamma_HA .- GR.*TL_HA; x=x+n;
    # TL_HB
    dydt[x:x+n-1] = beta_B_vec.*m_HB.*R .- ku_TL*TL_HB .- gamma_HB .- GR.*TL_HB; x=x+n;
    # TL_HC
    dydt[x:x+n-1] = beta_C_vec.*m_HC.*R .- ku_TL*TL_HC .- gamma_HC .- GR.*TL_HC; x=x+n;

    # R
    dydt[x:x+n-1] = gamma_R + ku_TL*(TL_all -(TL_HA + TL_HB + TL_HC))- kb_TL*(m_all - (m_HA + m_HB + m_HC)).*R + gamma_all
     + ku_TL*TL_HA - beta_A_vec.*m_HA.*R + ku_TL*TL_HB - beta_B_vec.*m_HB.*R + ku_TL*TL_HC - beta_C_vec.*m_HC.*R - GR.*R;
    x = x+n;

    # Z
    dydt[x:x+n-1] = gamma_Z - GR.*Z; x=x+n;
    # Q
    dydt[x:x+n-1] = gamma_Q - GR.*Q; x=x+n;
    
    # HA
    dydt[x:x+n-1] = gamma_HA - (GR.+p_deg).*HA; x=x+n;
    # HB
    dydt[x:x+n-1] = gamma_HB - (GR.+p_deg).*HB; x=x+n;
    # HC
    dydt[x:x+n-1] = gamma_HC - (GR.+p_deg).*HC;

    
end