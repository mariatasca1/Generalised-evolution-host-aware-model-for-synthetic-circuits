# Additional script that labels each state systematically with a
# coordinate ('coord') of length 's', then creates the cell array
# state_connections that collects info on how each state is connected to
# one another. This is required for automatically assigning 'part' and
# 'z' values to each state in the correct order. This script works for
# "doi.org/10.1101/2023.04.08.536106".



# Arguments
# `s::Int`: The size of each dimension.
# `d::Int`: The number of dimensions.
# `n::Int`: The total number of states.

# Returns
#`state_connections::Array{Any}`: A 2D array where each row represents a state and contains:
#    1. Index of the state.
#    2. Coordinate of the state.
#    3. List of coordinates of upstream states.
#    4. List of index pairs specifying mutation probabilities for upstream states.
#    5. List of coordinates of downstream states.
#    6. List of index pairs specifying mutation probabilities for downstream states.

function f_allocate_coords(s,d,n)

# Create systematic indices for each state

coord_system = zeros(n,d); # Initialise

for i = 1:d
    
    # Goal: create column of indices whose rows represent each state, one
    # dimension at a time. Order this top-down. E.g. for [s3,d2]: column
    # vectors would be:[2,1,0,2,1,0,2,1,0], [2,2,2,1,1,1,0,0,0].
    
    # 1. We want to cycle coords every s^(i-1)
    # -> x = repelem(s-1:-1:0, s^(i-1))
    # e.g. when s=3,d=3, cycle every 3^2=9
    
    # 2. We then repeat THAT cycle s^(d-i) times
    # repmat(x, [1,s^(d-i)])
    # e.g. for s=3,d=4: when d=2 we repeat s^(4-2)=9 times
    
    # 3. Order coords in reverse: a:-1:b
    
    # Break down the complex expression into intermediate steps
    cycle_length = s^(i-1);  # Length of each cycle
    repeat_count = s^(d-i) ; # Number of times to repeat the cycle
    base_sequence = s-1:-1:0 ; # Base sequence for coordinates
    repeated_sequence = repeat(base_sequence, inner =  cycle_length) ; # Create repeated sequence
    coord = repeat(repeated_sequence, outer =  repeat_count);  # Final coordinate array
    coord_system[:,i] = coord;
end



# Create 'state_connections' cell array

# Rows: states in order. Columns: (1) index of state. (2) coord of state.
# (3) List of coords of upstream states. (4) List of index pairs that
# each specify a value from the z_matrix for the mutation probability of
# an upstream state. (5-6) Same as (3-4) but for downstream states.

state_connections = Array{Any}(undef, n,6); # Initialise

# First two rows are simple
for i = 1:n
    state_connections[i,1] = i; # Column? 1
    state_connections[i,2] = coord_system[i,:]; # Column? 2
end

# 'From' (upstream) and 'To' (downstream coords are more involved.
# Here-in, 'partial coord' (pCoord) refers to a single digit within the
# full coord. CurrentState = s - pCoord

for i = 1:n # For every state
    state_connections[i, 3] = Vector{Any}()
    state_connections[i, 4] = Vector{Any}()
    state_connections[i, 5] = Vector{Any}()
    state_connections[i, 6] = Vector{Any}()
    for j = 1:d # Check transitions for each dimension
        
        old_pCoord = state_connections[i,2][j]; # Access the j-th element of the vector stored in the cell
        
        # 1. Find all UPSTREAM states        
        new_pCoord = old_pCoord; # Take coordinate of current subpop and record pCoord of current dim
        
        while new_pCoord != s-1 # Are backwards transitions possible?            
            new_pCoord = new_pCoord + 1; # For next subpop, increment state
            pre_subpop = copy(state_connections[i,2]); # Initialise coord for upstream subpop
            pre_subpop[j] = new_pCoord; # Set coord
            push!(state_connections[i,3], pre_subpop) # Add to cell array   
            
            # For upstream, z_type relates to *current* subpop
            z_type = [old_pCoord+1, j]; # Rows = decreasing severity. Index = pCoord+1
            push!(state_connections[i,4], z_type);
        end
        
        # 2. Find all DOWNSTREAM states
        new_pCoord = old_pCoord;
        
        while new_pCoord != 0            
            new_pCoord = new_pCoord - 1; # For next subpop, reduce state
            post_subpop = copy(state_connections[i,2]);
            post_subpop[j] = new_pCoord;
            push!(state_connections[i,5], post_subpop);
            
            # For downstream, z_type relates to *new* subpop
            z_type = [new_pCoord+1, j];
            push!(state_connections[i,6], z_type);
        end
        
    end # Dimensions
    
end # States

return state_connections
end