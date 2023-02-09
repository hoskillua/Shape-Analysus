using Makie
using Printf

include("utils.jl")

filename = "moomoo.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

M = getMeshData(X, T)

## Define distributions; match one vertex to the entire mesh
source = 1
q = zeros(M.numVertices, M.numVertices)
q[source,:] .= 1 ./ M.areaWeights[source]
p = diagm(0 => 1 ./ M.areaWeights)

## Regularized EMD code
tol = 1e-6
alpha = .00001
steps = 3
kernel(x) = heatDiffusion(x,M,alpha,steps)
kernelTranspose(x) = heatDiffusion(x,M,alpha,steps;transpose=true)

p = p .+ eps()
q = q .+ eps()

# You'll compute geodesic distances from point 1 to all other points
function emd(p, q, kernel, kernelTranspose)
    dists = zeros(size(p, 2), 1)

    niter = 100
    v = ones(size(p))
    w = ones(size(q))

    for i=1:niter
        ## YOUR CODE HERE TO UPDATE v, w, VARIABLES ##
        # Don't forget to multiply by area weights

        ## END CODING ASSIGNMENT ##

        oldDistances = dists

        ll(x) = real.(log.(Complex.(x)))
        ## YOUR CODE HERE TO COMPUTE DISTANCE USING v, w ##

        ## END CODING ASSIGNMENT

        println(size(dists))
        change = norm(oldDistances-dists,2)
        @printf("Iteration %d: %f\n", i, change)
        if change < tol
            break
        end
    end

    dists = sqrt.(max.(dists,0.0))
    return dists
end

dists = emd(p, q, kernel, kernelTranspose)
showdescriptor(X, T, dists)
