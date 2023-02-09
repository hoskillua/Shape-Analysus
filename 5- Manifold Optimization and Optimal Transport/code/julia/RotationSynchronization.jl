using LinearAlgebra
using Manifolds
using Manopt

function RotationSynchronization(Oij)

    N = Int(size(Oij, 1) / 3)

    ## Riemannian optimization
    M = PowerManifold(Rotations(3), N);

    function mycost(M, Y)
        ## Homework: implement the cost
        return 0
    end

    function mygrad(M, Y)
        ## Homework: implement the Riemannian gradient
    end

    function myhess(M, Y, V)
        ## Homework (optional): implement the Riemannian hessian-vector product
    end

    # Rhat = trust_regions(M, mycost, mygrad, myhess, random_point(M))
    Rhat = gradient_descent(M, mycost, mygrad, random_point(M))
    RijHat = reshape(Rhat, 3, 3 * N)' * reshape(Rhat, 3, 3 * N)

    return Rhat, RijHat

end
    








