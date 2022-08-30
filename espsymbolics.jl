using Symbolics
using Symbolics: scalarize

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

    # Check
    n > size(V, 1) && return zero(T)
    n == 0 && return one(T)
    size(V, 1) == 1 && return V[1]

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


"""
Keep ESP terms containing the input indices.

# Arguments
- `V`: vector of dimension m
- `n`: degree of ESP
- `vind`: V indices 
"""
function elemesp(V::AbstractVector{T},
                 n::Int, 
                 vind::Tuple) where T<:Num

    # Check
    size(vind, 1) > n && return zero(T)

    # Initialize
    m = size(V, 1)
    LV = zeros(T, m)
    
    # Modify
    for i = 1:m
        LV[i] = V[i] 
    end
    for p = 1:size(vind, 1)
        LV[vind[p]] = zero(T) 
    end

    # Compute
    indval = scalarize(sumesp(LV, n - size(vind, 1)))
    for p = 1:size(vind, 1)
	indval *= V[vind[p]] 
    end

    return indval 
end


let
    
    # Print
    println("Symbolic ESP")

    # User inputs
    m = 8
    n = 4
    p = (4, 3, 5)

    # Variables
    @variables V[1:m] 

    # Get ESP
    Emn = expand(sumesp(V, n))

    # Element extract
    Emnp = expand(elemesp(V, n, p))

    # Print
    println("---------------------")
    println("m: $m, n: $n")
    println("---------------------")
    println("E$m$n: $Emn")
    println("---------------------")
    println("E$m$n -> indices $p: $Emnp")
    println()

    # Exit	
    nothing	
end	


