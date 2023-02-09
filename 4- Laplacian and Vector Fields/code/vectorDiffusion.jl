using Arpack
using SparseArrays
using Statistics
using Rotations
using WGLMakie

include("utils/utils.jl")

function vectorDiffusion()
    meshfile = "meshes/moomoo.off"
    X, T = Utils.readoff(meshfile)
    nv = size(X, 1)
    nt = size(T, 1)

    data = Utils.MeshData(X,T);

    ## PROBLEM 3(a)
    ### YOUR CODE HERE - implement the connection laplacian operator
    # data: mesh data structure
    # CL: |3V| x |3V| "sparse connection laplacian operator"
    function buildConnectionLaplacian(data)
        CL = spzeros(3 * data.nv, 3 * data.nv)
        return CL
    end
    ### END HOMEWORK PROBLEM

    CL = buildConnectionLaplacian(data)

    # build vector valued signal: delta vector heat source at highest Y axis point.
    dt = sqrt(mean(data.edgeLengths));
    signal = zeros(3, data.nv)
    maxind = argmax(data.verts[:, 2])
    t1 = cross(data.vertNormals[maxind, :], randn(3))
    t1 = t1 / norm(t1)
    signal[:, maxind] = t1
    signal = signal[:]

    ### PROBLEM 3(b) - VECTOR HEAT METHOD
    # todo: compute diffused signal for various dt and normalize
    phi = zeros(data.nv, 3)
    ### END HOMEWORK PROBLEM

    ## Visualization
    fig, ax, msh = mesh(X, T, color=:black, shading=true, figure=(resolution=(1000,1000),))
    arrows!(Point3.(eachrow(X)), Point3.(eachrow(phi)), linewidth=200, linecolor=:red, arrowcolor=:red, arrowsize=0.2)#, lengthscale=0.1)
    arrows!([Point3(data.verts[maxind, :])], [Point3(t1)], linewidth=400, linecolor=:cyan, arrowcolor=:cyan, arrowsize=0.4, lengthscale=2)
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)
end

