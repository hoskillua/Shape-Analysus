# using Arpack
# using SparseArrays
# using Statistics
# using WGLMakie
using Manopt, Manifolds, Random

function karchergrad(M, x, Z)
    ### Homework: Implement the Riemannian gradient as rgrad
    rgrad = 0*x
    return rgrad
end

function karcherhess(M, x, v, Z)
    ### Homework: Implement the Riemannian hessian as rhess
    rhess = 0

    # Multiply Riemannian Hessian with v
    rhessv = rhess*v
    return rhessv
end

function SphereKarcherMeans()
    N = 5
    Z = abs.(randn(N,2)); Z = Z./sqrt.(sum(Z.*Z,dims=2));
    
    kcost(M, y) = sum(acos.(Z*y).^2)
    kgrad(M, x) = karchergrad(M, x, Z)
    khess(M, x, v) = karcherhess(M, x, v, Z)

    M = Manifolds.Sphere(1)
    x0 = [-1, -1]/sqrt(2)

    xstar2 = gradient_descent(M, kcost, kgrad, x0)
    # xstar1 = trust_regions(M, kcost, kgrad, khess, x0)
end


    








