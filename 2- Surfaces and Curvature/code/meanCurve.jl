using Arpack
using Printf
using SparseArrays

include("utils.jl")

filename = "meshes/moomoo.off"
# filename = "meshes/166.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

# ADD CODE TO COMPUTE SURFACE AREA HERE #############
function surfaceArea(X, T)
    return 0
end
# END HOMEWORK ASSIGNMENT #########

@printf("The surface area of %s is %f\n", filename, surfaceArea(X, T))

# ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE ##########
function cotLaplacian(X, T)
    return spzeros(nv, nv)
end
# END HOMEWORK ASSIGNMENT #######

# Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X, T)
eigenvals, _ = eigs(L, nev=10, which=:SM)
println(eigenvals)
println(norm(L - L'))

# ADD CODE FOR DIVIDED DIFFERENCES HERE ######
function dividedDifferences(X, T)
    gradApprox = zeros(nv, 3);
    return gradApprox
end
# END HOMEWORK ASSIGNMENT #########

# Check that gradApprox and .5*L*X are similar
println(norm(.5 .* (L*X) - dividedDifferences(X,T)))

# ADD CODE FOR COMPUTING THE BARYCENTRIC AREA VECTOR HERE ######
function barycentricArea(X, T)
    return zeros(nv)
end
# END HOMEWORK ASSIGNMENT ########

# ADD CODE FOR COMPUTING POINTWISE MEAN CURVATURE HERE #####
function meanCurvature(X, T)
    H = rand(nv)
    return H
end
# END HOMEWORK ASSIGNMENT #######

H = meanCurvature(X, T)
scene = showdescriptor(X, T, H)


## Mean curvature flow ##
function curvatureFlowEuler(X, T)
    Xt = copy(X)
    maxiters = 1000
    # ADD CODE FOR THE EXPLICIT INTEGRATOR HERE ####
    for t=1:maxiters
        Xt[:] = copy(Xt)
    end
    # END HOMEWORK ASSIGNMENT #####
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    # scene = showdescriptor(X, T, H)
end

function curvatureFlowImplicit(X, T)
    Xt = copy(X)
    maxiters = 1000
    # ADD CODE FOR SEMI-IMPLICIT INTEGRATOR HERE ####
    for t=1:maxiters
        Xt[:] = copy(Xt)
    end
    # END HOMEWORK ASSIGNMENT
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    # scene = showdescriptor(X, T)
end

# Uncomment either the implicit or explicit flow to see your results
# curvatureFlowEuler(X,T)
# curvatureFlowImplicit(X,T)
