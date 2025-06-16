# Additional script that converts information from 'state_connections'
# into more useable variables. The ODE for each state involves terms
# regarding (i) number of cells in connected upstream states, (ii)
# mutation probabilities from connected upstream states, (iii) mutation
# probabilities to downstream states. The cell array 'state_connections'
# has this info in an over-extensive and coded format. Hence,
# f_extract_state_connections transfers the required info into defined
# variables and converts 'coordinates' to numbered indices and
# 'z_matrix index pairs' into probabilities. For more contextual info,
# see "doi.org/10.1101/2023.04.08.536106".
# Code author: Duncan Ingram, 2023.

function f_extract_state_connections(state_connections, active_z_matrix, s, d, n)
      
# Pre-allocate variables
states_up = [Int[] for _ in 1:n];
z_values_up = [Float64[] for _ in 1:n];
states_dn = [Int[] for _ in 1:n];
z_values_dn = [Float64[] for _ in 1:n];
include("f_coord_converter.jl")

# Find upstream/downstream. For each state, store indicies of their up/dn
# states and z-values. Info for 'cells_dn' is not needed, but can be
# stored regardless.

for i=1:n

    num_up = size(state_connections[i,3], 1);
    num_dn = size(state_connections[i,5], 1);
    
    # Upstream
    for j = 1:num_up
        
        # States
        coord_up = vec(state_connections[i,3][j]); # Get coord of first upstream state
        idx = f_coord_converter(coord_up, s, d); # What idx is that state?
        push!(states_up[i], idx);
        
        # z: values from state_connections act as indicies for z_mat
        state = state_connections[i,4][j][1];
        dim = state_connections[i,4][j][2];        
        push!(z_values_up[i], active_z_matrix[Int(state), Int(dim)]);
    end    
    
    # Downstream
    for j = 1:num_dn
        
        # " "
        coord_dn = vec(state_connections[i,5][j]);
        idx = f_coord_converter(coord_dn, s, d);
        push!(states_dn[i], idx);
        
        # " "
        state = state_connections[i,6][j][1];
        dim = state_connections[i,6][j][2];        
        push!(z_values_dn[i], active_z_matrix[Int(state),Int(dim)]);
    end    
end

# Return a list containing:
# - states_up: indices of upstream states for each state
# - z_values_up: z-matrix values corresponding to upstream states
# - z_values_dn: z-matrix values corresponding to downstream states
return [states_up, z_values_up, z_values_dn]
end