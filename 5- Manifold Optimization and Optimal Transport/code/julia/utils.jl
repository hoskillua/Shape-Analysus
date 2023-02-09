using SparseArrays
using LinearAlgebra
using Makie

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

"""
    norms(X; dim)
"""
function norms(X::Matrix{Float64}; dim=1)
    if dim == 1
        return [norm(X[i,:]) for i=1:size(X,1)]
    else
        return [norm(X[:,i]) for i=1:size(X,2)]
    end
end

"""
    cotLaplacian(X, T)
"""
function cotLaplacian(X::Matrix{Float64},
                      T::Matrix{Int})
    nv = size(X,1)
    nt = size(T,1)

    I = T[:,1]; J = T[:,2]; K = T[:,3];
    crossprods = zeros(nt,3)
    for t=1:nt
        crossprods[t,:] = cross(X[J[t],:]-X[I[t],:], X[K[t],:]-X[I[t],:])
    end
    areas = .5*norms(crossprods,dim=1)

    cots = zeros(nt,3)
    for i=0:2
        I = T[:,i+1]; J = T[:,(i+1)%3+1]; K = T[:,(i+2)%3+1];
        e1 = X[J,:]-X[I,:]
        e2 = X[K,:]-X[I,:]
        cots[:,i+1] = sum(e1.*e2,dims=2)./(2*areas)
    end

    I = T[:,1]; J = T[:,2]; K = T[:,3];
    L = sparse([I;J;K], [J;K;I], [cots[:,3];cots[:,1];cots[:,2]], nv, nv)
    L = L + L'
    rowsums = vec(sum(L,dims=2))
    L = L - spdiagm(0 => rowsums)
    return -.5 * L
end

"""
    barycentricAreas(X,T)
"""
function barycentricAreas(X::Matrix{Float64},
                          T::Matrix{Int64})
    nv = size(X,1)
    nt = size(T,1)

    I = T[:,1]; J = T[:,2]; K = T[:,3];
    crossprods = zeros(nt,3)
    for t=1:nt
        crossprods[t,:] = cross(X[J[t],:]-X[I[t],:], X[K[t],:]-X[I[t],:])
    end
    areas = .5*norms(crossprods,dim=1)

    baryAreas = zeros(nv)
    for t=1:nt
        baryAreas[T[t,1]] += areas[t]/3
        baryAreas[T[t,2]] += areas[t]/3
        baryAreas[T[t,3]] += areas[t]/3
    end

    return baryAreas
end

"""
    massMatrix(X,T)
"""
function massMatrix(X::Matrix{Float64},
                    T::Matrix{Int64})
    return spdiagm(0 => barycentricAreas(X,T))
end

mutable struct Mesh
    cotLaplacian
    areaWeights
    numVertices::Int64
    numTriangles::Int64
    vertices::Matrix{Float64}
    triangles::Matrix{Int64}
end

"""
    getMeshData(X,T)
"""
function getMeshData(X::Matrix{Float64},
                     T::Matrix{Int64})
    nv = size(X, 1)
    nt = size(T, 1)
    L = cotLaplacian(X, T)
    A = barycentricAreas(X, T)

    mesh = Mesh(L, A, nv, nt, X, T)
    return mesh
end

"""
    heatDiffusion(signal, mesh, time, steps)
"""
function heatDiffusion(signal, mesh, time, steps; transpose=false)
    h = time / steps
    nv = mesh.numVertices

    blurInverse = spdiagm(0 => mesh.areaWeights) - h * mesh.cotLaplacian
    result = copy(signal)
    if transpose
        for i = 1:steps
            result = (blurInverse' \ result) .* mesh.areaWeights
        end
    else
        for i = 1:steps
            result = blurInverse \ (result .* mesh.areaWeights)
        end
    end

    return result
end
