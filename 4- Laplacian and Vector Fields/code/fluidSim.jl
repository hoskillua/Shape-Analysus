using Arpack
using SparseArrays
using SuiteSparse
using Statistics
using BenchmarkTools
using WGLMakie

include("utils/utils.jl")
include("getDivGrad.jl")

function fluidSim()
    meshfile = "meshes/spherers.off"
    X, T = Utils.readoff(meshfile)
    nv = size(X, 1)
    nt = size(T, 1)

    data = Utils.MeshData(X,T);
    Div, Grad = getDivGrad(data)
    meanEdgeLength = mean(data.edgeLengths)
    dt = 0.0005

    # initialize velocity field and timestep
    V = [zeros(data.nf) -data.faceCenters[:, 3] data.faceCenters[:, 2]]
    V = V - sum(V .* data.faceNormals, dims=2) .* data.faceNormals
    V = V ./ norm.(eachrow(V))
    V .= (10 * meanEdgeLength) .* V
    div = zeros(data.nv)
    pressure = zeros(data.nv)

    # initialize color field
    color = cos.(5 .* data.verts[:, 3])

    # initialize source region
    mx = maximum(data.faceCenters[:, 3])
    sourceInd = findall(mx .- data.faceCenters[:, 3] .< 0.05)


    # Visualize starting conditions
    Vvis = Node(Point3.(eachrow(V * (2 * meanEdgeLength / maximum(norm.(eachrow(V)))))))
    cvis = Node(color)
    fig, ax, msh = mesh(X, T, color=cvis, shading=false, figure=(resolution=(1000,1000),))
    arrows!(Point3.(eachrow(data.faceCenters)), Vvis, linewidth=2, linecolor=:red, arrowcolor=:red, arrowsize=0.005, lengthscale=0.5)
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)

    # Factorize laplacian
    decL = ldlt(Symmetric(data.cotLaplacian))

    # fluid simulation loop.
    maxiters = 100000;
    t = 0
    for i=1:maxiters
        sourceForce = (0.2 * meanEdgeLength / dt) * [cos(0.2 * pi * t) sin(0.2 * pi * t) 0]
        sourceForce = sourceForce .- sum(sourceForce .* data.faceNormals, dims=2) .* data.faceNormals

        ### PROBLEM 4(d) - YOUR CODE HERE
        # integrate source force
        
        ### PROBLEM 4(e) - YOUR CODE HERE
        # Implement pressure projection. Make sure V is divergence free
        
        # Update visualization
        if mod(i, 66) == 0
            Vvis[] = Point3.(eachrow(V * (2 * meanEdgeLength / maximum(norm.(eachrow(V))))))
            cvis[] = color
            sleep(0.01)
        end
        
        ### PROBLEM 4(a) - YOUR CODE HERE
        # construct Dv operator

        ### PROBLEM 4(b) - YOUR CODE HERE
        # advect color field

        ### PROBLEM 4(c) - YOUR CODE HERE
        # implement fluid self-advection

        ### END HOMEWORK PROBLEM
        
        t += dt
    end
end


