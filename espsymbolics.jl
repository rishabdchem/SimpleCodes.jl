using Symbolics

#=====================================================
Symbolical evaluation of an m-variable
n-degree elementary symmetric polynomial.
The integers m and n are user defined.
=====================================================#

"""
SumESP algorithm.
See DOI: 10.1016/j.amc.2015.08.134.

# Arguments
- `V`: vector of dimension m
- `n`: degree of ESP
"""
function sumesp(V::AbstractVector{T}, 
	        n::Int) where T<:Num

    # Initialize
    m = size(V, 1)
    E = zeros(T, m, n+1)
    E[:, 1] = ones(T, m)
    E[1, 2] = V[1]
    
    # SumESP
    for i = 2:m, j = max(1, i + n - m):min(i, n)
        E[i, j+1] = E[i-1, j+1] + V[i] * E[i-1, j]
    end
    
    return E[m, n+1] 
end


let
    
    # Print
    println("Symbolic ESP")

    # User inputs
    m = 8
    n = 4

    # Variables
    @variables V[1:m] 

    # Get ESP
    Emn = expand(sumesp(V, n))

    # Print
    println("m: $m, n: $n")
    println("---------------------")
    println("ESP is: $Emn")
    println()

    # Exit	
    nothing	
end	


