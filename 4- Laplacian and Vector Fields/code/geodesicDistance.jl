using Arpack
using SparseArrays
using Statistics
using WGLMakie

include("utils/utils.jl")
include("getDivGrad.jl")

function geodesicDistance()
    meshfile = "meshes/moomoo.off"
    X, T = Utils.readoff(meshfile)
    nv = size(X, 1)
    nt = size(T, 1)

    # load mesh
    data = Utils.MeshData(X,T);

    # build signal: delta heat source at highest Y axis point. This defines the source vertex.
    dt = sqrt(mean(data.edgeLengths));
    signal = zeros(nv,1);
    maxind = argmax(data.verts[:, 2]);
    signal[maxind] = 1;
    ## The heat method 

    ### PROBLEM 2 - THE HEAT METHOD
    ### YOUR CODE HERE - short-time heat diffusion
    phi = zeros(data.nv)

    Div, Grad = getDivGrad(data)

    # YOUR CODE HERE - compute the normalized gradient of phi
    qvec = zeros(data.nf, 3)

    # YOUR CODE HERE - compute divergence of normalized gradient
    divVals = zeros(data.nv)

    # YOUR CODE HERE - final geodesic distances
    dist = zeros(data.nv)
    ### END HOMEWORK PROBLEM

    ## Visualization
    fig, ax, msh = mesh(X, T, color=phi[:], shading=true, figure=(resolution=(1000,500),))
    arrows!(Point3.(eachrow(data.faceCenters)), Point3.(eachrow(qvec)), linewidth=200, linecolor=:red, arrowcolor=:red, arrowsize=0.1)#, lengthscale=0.1)
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax

    ax2, msh2 = mesh(fig[1, 2], X, T, color=divVals[:], shading=false)
    oax2 = ax2.scene[OldAxis]
    oax2.showgrid = (false, false, false)
    oax2.showticks = (false, false, false)

    ax3, msh3 = mesh(fig[1, 3], X, T, color=dist[:], shading=false, colormap=:prism)
    oax3 = ax3.scene[OldAxis]
    oax3.showgrid = (false, false, false)
    oax3.showticks = (false, false, false)

    display(fig)
end
