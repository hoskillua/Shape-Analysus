module Utils

using LinearAlgebra
using SparseArrays
using Statistics
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

mutable struct MeshData
    verts :: Matrix{Float64}
    faces :: Matrix{Int}
    edges :: Matrix{Int}
    face2edge :: Matrix{Int}
    # faceOrientation :: Array{Int, 2}
    nv :: Int
    ne :: Int
    nf :: Int
    faceAreas :: Vector{Float64}
    faceNormals :: Matrix{Float64}
    faceCenters :: Matrix{Float64}
    vertAreas :: Vector{Float64}
    vertNormals :: Matrix{Float64}
    edgeLengths :: Vector{Float64}
    primalIncidence :: SparseMatrixCSC{Float64, Int}
    cotweights :: SparseMatrixCSC{Float64, Int}
    cotLaplacian :: SparseMatrixCSC{Float64, Int}
    massMatrix :: SparseMatrixCSC{Float64, Int}
    FtoV :: SparseMatrixCSC{Float64, Int}
    VtoF :: SparseMatrixCSC{Float64, Int}

    @views function MeshData(X, T)
        mesh = new(X, T)
        mesh.nv = size(X, 1)
        mesh.nf = size(T, 1)

        areaNormals = reduce(hcat, cross.(eachrow(X[T[:, 2], :] - X[T[:, 1], :]), eachrow(X[T[:, 3], :] - X[T[:, 1], :])))'
        mesh.faceAreas = 0.5 .* norm.(eachrow(areaNormals))
        mesh.faceNormals = areaNormals ./ (2 .* reshape(mesh.faceAreas, :, 1))
        mesh.vertAreas = zeros(mesh.nv)
        for f = 1:mesh.nf
            for i = 1:3
                mesh.vertAreas[T[f, i]] += mesh.faceAreas[f] / 3
            end
        end
        mesh.faceCenters = reshape(mean(X[T, :], dims=2), mesh.nf, 3)

        faceEdges = sort(reshape(T[:, [1 2; 2 3; 3 1]], :, 2), dims = 2)
        mesh.edges = unique(faceEdges, dims=1)
        mesh.ne = size(mesh.edges, 1)
        mesh.face2edge = reshape(indexin(eachrow(faceEdges), collect(eachrow(mesh.edges))), mesh.nf, 3)
        # mesh.faceOrientation = sign(mesh.edges[mesh.face2edge, 2] - mesh.edges[mesh.face2edge, 1])
        mesh.edgeLengths = norm.(eachrow(X[mesh.edges[:, 2], :] - X[mesh.edges[:, 1], :]))
        
        mesh.primalIncidence = sparse(vec(repeat(1:mesh.ne, 1, 2)), vec(mesh.edges), vec(repeat([-1 1], mesh.ne, 1)), mesh.ne, mesh.nv)

        eij = reshape(X[T[:, [2, 3, 1]], :] - X[T, :], :, 3)
        eik = reshape(X[T[:, [3, 1, 2]], :] - X[T, :], :, 3)
        cotangents = reshape(dot.(eachrow(eij), eachrow(eik)) ./ norm.(cross.(eachrow(eij), eachrow(eik))), :, 3)
        mesh.cotweights = sparse(vec(mesh.face2edge[:, [2 3 1]]), vec(mesh.face2edge[:, [2 3 1]]), vec(cotangents), mesh.ne, mesh.ne)

        mesh.cotLaplacian = 0.5 * (mesh.primalIncidence' * mesh.cotweights * mesh.primalIncidence)
        mesh.massMatrix = Diagonal(mesh.vertAreas)
        mesh.FtoV = mesh.massMatrix \ sparse(vec(T), vec(repeat(1:mesh.nf, 1, 3)), vec(repeat(mesh.faceAreas ./ 3, 1, 3)), mesh.nv, mesh.nf)
        mesh.VtoF = sparse(vec(repeat(1:mesh.nf, 1, 3)), vec(T), fill(1/3, length(T)), mesh.nf, mesh.nv)
        mesh.vertNormals = mesh.FtoV * mesh.faceNormals
        mesh.vertNormals .= mesh.vertNormals ./ sqrt.(sum(mesh.vertNormals.^2, dims=2))

        return mesh
    end
end

end