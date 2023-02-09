clear; clc;
addpath('utils/');

%% Load mesh
filename = 'moomoo.off';
[X,T] = readOff(filename);
nv = size(X,1);
nt = size(T,1);
mesh = getMeshData(X,T,10);

%% Define distributions
source = 1;
q = zeros(mesh.numVertices,mesh.numVertices);
q(source,:) = 1./mesh.areaWeights(source);
p = spdiags(1./mesh.areaWeights,0,mesh.numVertices,mesh.numVertices);

%% Regularized EMD code
tol = 1e-6;
alpha = .00001;
steps = 3;
kernel = @(x) heatDiffusion(x,mesh,alpha,steps);
kernelTranspose = @(x) heatDiffusion(x,mesh,alpha,steps,1);

p = p+eps;
q = q+eps;

% You'll compute geodesic distances from point 1 to all others
distances = zeros(size(p,2),1);

niter = 30;
v = ones(size(p));
w = ones(size(q));

aw = bsxfun(@times,mesh.areaWeights,w);
for i=1:niter
    % YOUR CODE HERE TO UPDATE v, w VARIABLES %%%
    % Remember to multiply by areaWeights

    % END CODING ASSIGNMENT %%%
    
    oldDistances = distances;

    ll = @(x) real(log(x));
    % YOUR CODE HERE TO COMPUTE DISTANCE USING v, w %%%
    
    % END CODING ASSIGNMENT %%%
    
    change = norm(oldDistances-distances,'fro');
    fprintf('Iteration %d: %g\n', i, change);
    if change<tol
        break;
    end
end

distances = sqrt(max(distances,0));
showDescriptor(mesh,distances);