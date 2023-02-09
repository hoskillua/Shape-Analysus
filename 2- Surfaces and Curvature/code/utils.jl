using LinearAlgebra
using WGLMakie

"""
    readoff(filename)

Read a .off file and return a list of vertex positions and a triangle matrix.
"""
function readoff(filename::String)
    X, T = open(filename) do f
        s = readline(f)
        if s[1:3] != "OFF"
            error("The file is not a valid OFF one.")
        end
        s = readline(f)
        nv, nt = parse.(Int, split(s));
        X = zeros(Float64, nv, 3);
        T = zeros(Int, nt, 3);
        for i=1:nv
            s = readline(f)
            X[i,:] = parse.(Float64, split(s))
        end
        for i=1:nt
            s = readline(f)
            T[i,:] = parse.(Int64, split(s))[2:end] .+ 1
        end
        X, T
    end
end

"""
    showmesh(X, T)

Plot a mesh represented by coordinates `X` and connectivity `T` and colored by z value.
If `frame` is set to true, also draw a wireframe of the mesh.

See also: [`showdescriptor`](@ref)
"""
function showmesh(X::Matrix{Float64}, T::Matrix{Int};
                  frame=false)
    Xmut = Node(X)
    ymut = lift(X -> @view(X[:, 2]), Xmut)
    fig, ax, msh = mesh(Xmut, T, color=ymut, shading=true, figure=(resolution=(1000, 1000),))
    if frame
        wireframe!(msh[1], color=(:black, 0.6), linewidth=50)
    end
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)
    Xmut
end

"""
    updatemesh!(X, T)

Update the current display window with new mesh. Use this for when
you need to visualize an animation.
"""
function updatemesh!(Xmut::Observable, X::Matrix{Float64})
    Xmut[] = X
end

"""
    showdescriptor(X, T, desc)

Plot a mesh represented by coordinates `X` and connectivity `T` and colored by the
value of `desc`. If `frame` is set to true, also draw a wireframe of the mesh.
"""
function showdescriptor(X::Matrix{Float64},
                        T::Matrix{Int},
                        desc::Vector{Float64};
                        frame=false)
    if length(desc) != size(X, 1)
        error("There must be a descriptor for each vertex in the mesh")
    end
    
    fig, ax, msh = mesh(X, T, color=desc, shading=true, figure=(resolution=(1000, 1000),))
    if frame
        wireframe!(msh[1], color=(:black, 0.6), linewidth=50)
    end
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)
end
