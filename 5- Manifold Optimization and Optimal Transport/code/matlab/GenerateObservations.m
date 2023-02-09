function [Oij, RijTrue] = GenerateObservations(N, sigma)

%% Generate true rotations
Ri = randrot(3, N);
Ri = reshape(Ri, 3, 3 * N);
RijTrue = Ri' * Ri;

%% Add random noise
Oij = RijTrue + sigma * randn(size(RijTrue));
Oij = triu(Oij, 1) + triu(Oij, 1)' + diag(diag(Oij));

end