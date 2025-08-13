# Additional script that obtains the 'active' parameters for part
# strengths and z_values, and distributes them across the states in the
# correct order. 


#Constructs active parameters for part strengths and mutation probabilities 
#based on the active dimensions and distributes them across states in the 
#correct order.

# Arguments
# `active_dims::Vector{Int}`: A vector indicating active dimensions (1 for active, 0 for inactive).
# `s::Int`: Number of states.
# `d::Int`: Number of dimensions.
# `n::Int`: Number of constructs.
# `part_matrix::Matrix{Float64}`: Matrix containing part strengths.
# `z_matrix::Matrix{Float64}`: Matrix containing mutation probabilities.

# Returns
# `active_z_matrix::Matrix{Float64}`: Matrix of active mutation probabilities.
# matrix::Matrix{Float64}`: Matrix of active part strengths distributed across constructs.

function f_active_construct_params(active_dims, s, d, n, part_matrix, z_matrix)
                                       
# 1. Create arrays of values associated with the active dimensions

# Part strengths
active_part_matrix = zeros(s,d);
j = length(active_dims)
count = 1;
for i = 1:j
    if active_dims[i] == 1
        if s==1 # Different indexing if all E cells
            active_part_matrix[1,count] = part_matrix[1,i];
        else
            @show part_matrix[:,i]
            active_part_matrix[1:s,count] = part_matrix[vcat(1:s-1, size(part_matrix,1)), i];
            # E.g. if s=3, parts_mat_all indicies are 1,2,4
        end
        count=count+1;
    end
end

# Mutation probabilities
active_z_matrix = zeros(s-1,d);
count = 1;
for i = 1:j
    if active_dims[i] == 1
        if s == 2
            active_z_matrix[1,count] = z_matrix[1,i];
            count=count+1;
        else
            active_z_matrix[1:s-1,count] = z_matrix[(1:s-2,size(z_matrix, 1)), i];
            # E.g. if s=3, z_mat_all indicies are 1,3
            count=count+1;
        end
    end
end


# 2. Allocate active parts (cols) to subpops in order (rows)

parts_combos = zeros(n,d);

for i = 1:d
    # For explanation, see f_allocate_coords.m
    @show active_part_matrix[1:end, i]
    repeated_elements = repeat(active_part_matrix[1:end, i], inner = s^(i - 1));
    @show repeated_elements
    @show s^(d-i)
    parts_combos[:, i] = repeat(repeated_elements, outer = s^(d - i));
end


# 3. Convert each column into correct part variable

# This section assigns the correct values to each active part manually
# and in order of the 'active_dims' vector [gene_A_prom, gene_A_RBS,
# gene_B_prom etc.]. It has code to account for up to three seperate
# genes. For >3 genes, more code is required in the same format as below.
matrix = zeros(Float64, length(active_dims), n)
ActiveDim = 1; # For allocating active parts
x = 1; # Part Type no.
num_genes = length(active_dims)/2; # '2': we consider prom/RBS for each gene
#each gene will have promoter and RBS so the length of the vector will always be even

############ Constructs ##########
for i in 1: Int(num_genes)
    if active_dims[x] == 1 #Promoter
        
        alpha_vec = transpose(parts_combos[:,ActiveDim]);
        ActiveDim = ActiveDim  +1;
    else
        alpha_vec = fill(part_matrix[1,x],1 , n);
    end
    
    matrix[2*i-1, :] = alpha_vec;
    x = x+1;

    if active_dims[x] == 1 # RBS
        beta_vec = transpose(parts_combos[:, ActiveDim]);
        ActiveDim = ActiveDim + 1;
    else
        beta_vec = fill(part_matrix[1,x], 1, n);
    end
    matrix[2*i, :] = beta_vec;
    x = x + 1;
end
    # `active_z_matrix`: A matrix containing active mutation probabilities for the constructs.
    # `matrix`: A matrix where each row represents a construct and columns represent active part strengths 
    #           distributed across constructs in the correct order.
    return active_z_matrix, matrix
end