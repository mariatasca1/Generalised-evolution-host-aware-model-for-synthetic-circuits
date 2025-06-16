# Additional script that converts a coordinate into an index. E.g. in a
# [s2,d3] framework, the coordinate [1,1,1] would translate to index '1',
# (i.e. the first state), and [0,0,0] to index '8' (i.e. the last state).
# For more contextual info, see "doi.org/10.1101/2023.04.08.536106".
#Code author: Duncan Ingram, 2023.

#Convert a coordinate into an index based on the [s2, d3] framework.

# Arguments
#`coord::Vector{Int}`: A vector representing the coordinate.
#`s::Int`: The base value for the coordinate system.
#`d::Int`: The dimensionality of the coordinate system.

# Returns
# `index::Int`: The calculated index corresponding to the given coordinate.

function f_coord_converter(coord, s, d)

    index = 1;
    powers_of_s = [s^(i-1) for i in 1:d]
    
    for i = 1:d
        unit_value = powers_of_s[i] * (s-coord[i]-1);
        index = index + unit_value;
    end

    return index
end