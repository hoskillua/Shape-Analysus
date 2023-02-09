include("utils/utils.jl")
using LinearAlgebra
using SparseArrays

## PROBLEM 2 - YOUR CODE HERE TO COMPUTE DIVERGENCE AND GRADIENT
# Div (|V| x 3|F| sparse) maps face-based vector fields to their vertex-based divergence
# Grad (3|F| x |V| sparse) maps a scalar function to its face-based gradient
function getDivGrad(data :: Utils.MeshData) :: Tuple{SparseMatrixCSC{Float64, Int}, SparseMatrixCSC{Float64, Int}}
    Div = spzeros(data.nv, 3 * data.nf)
    Grad = spzeros(3 * data.nf, data.nv)
    return Div, Grad
end
### END HOMEWORK PROBLEM