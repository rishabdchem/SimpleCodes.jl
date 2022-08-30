using Optim 
using Optim: only_fg!, minimum, minimizer

#=====================================================
General usage of the Optim.jl package.
=====================================================#

"""
Minimize a simple function with input parameters.

# Arguments
- `a`, `b`: parameters 
"""
function minimize(a::Float64,
	          b::Float64)
 
    # Guess	
    xguess = randn(2) 

    # Optimize with L-BFGS
    res = optimize(only_fg!( (f, G, x) -> getfg!(f, G, x, a, b) ), xguess, LBFGS())

    # Print
    println("---------------------")
    display(res)
    println()

    return minimum(res), minimizer(res) 
end


function getfg!(f, G, V, a, b)

    # Variables	
    x, y = V[1], V[2]

    # Gradient 
    if G !== nothing
        G[1] = 2 * x 
        G[2] = 2 * y
    end	
  
    # Function 
    if f !== nothing
        return (x ^ 2) + (y ^ 2) + a + b 
    end	
end


let
    
    # Print
    println("Minimization of a simple function")

    # User inputs
    a = rand(1.0:100.0) 
    b = rand(1.0:100.0) 

    # Minimize  
    f, V = minimize(a, b)

    # Print
    println("---------------------")
    println("a: $a, b: $b")
    println("Minimum: $f")
    println("Minimizer: $V")
    println()

    # Exit	
    nothing	
end	


