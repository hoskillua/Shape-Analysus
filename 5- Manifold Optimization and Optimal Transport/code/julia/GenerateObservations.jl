using LinearAlgebra
using Manifolds
using Manopt

function GenerateObservations(N :: Int, sigma :: Float64)
    # Generate true rotations
    Ri = random_point(PowerManifold(Rotations(3), N))
    Ri = reshape(Ri, 3, 3 * N)
    RijTrue = Ri' * Ri

    # Add random noise
    Oij = RijTrue + sigma * randn(size(RijTrue));
    Oij = triu(Oij, 1) + triu(Oij, 1)' + Diagonal(diag(Oij));

    return Oij, RijTrue
end