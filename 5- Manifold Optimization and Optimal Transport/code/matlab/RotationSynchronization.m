function [Rhat, RijHat] = RotationSynchronization(Oij)

N = size(Oij, 1) / 3;

%% Riemannian optimization
M = rotationsfactory(3, N);
problem.M = M;
problem.cost = @mycost;
problem.egrad = @myegrad;
problem.ehess = @myehess;

% figure; checkgradient(problem); pause;
% figure; checkhessian(problem); pause;

Rhat = trustregions(problem);

RijHat = reshape(Rhat, 3, 3 * N)' * reshape(Rhat, 3, 3 * N);

function [cost, store] = mycost(Y, store)
    %% Homework: implement the cost
    cost = 0;
end

function [egrad, store] = myegrad(Y, store)
    %% Homework: implement the Euclidean gradient as egrad
    egrad = zeros(size(Y));
end

function [ehessv, store] = myehess(Y, V, store)
    %% Homework: implement the Euclidean hessian-vector product as ehessv
    ehessv = zeros(size(V));
end

end
    








