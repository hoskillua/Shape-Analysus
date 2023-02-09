addpath('utils/');

meshfile = 'meshes/166.off';
[X,T] = readOff(meshfile);
nv = size(X,1);
nt = size(T,1);

A = surfaceArea(X,T);
fprintf('The surface area of %s is %f\n', meshfile, A);

% Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X,T);
[~,p] = chol(L);
fprintf('\nIf %d is 0, then L is PSD\n', p);
fprintf('Symmetry: %d\n', norm(L-L','fro'));

% Divided differences approximation
G = gradApprox(X,T);
fprintf('Difference between gradient and cotan Laplacian\n', norm(.5*L*X-G,'fro'));

H = meanCurvature(X,T);
showDescriptor(X,T,H);

%% Mean curvature flow
maxiters=1000;
Xt = X;
% ADD CODE FOR EXPLICIT INTEGRATOR HERE %%%%
for t=1:maxiters
    Xt(:) = X;
end
% END HOMEWORK ASSIGNMENT %%%%
% Uncomment to display mesh at the end
H = meanCurvature(Xt,T);
% showDescriptor(Xt,T,H);

Xt = X;
% ADD CODE FOR IMPLICIT INTEGRATOR HERE %%%%
for t=1:maxiters
    Xt(:) = X;
end
% END HOMEWORK ASSIGNMENT %%%%
% Uncomment to display mesh at the end
H = meanCurvature(Xt,T);
% showDescriptor(Xt,T,H);


%% Function definitions
% ADD CODE TO COMPUTE SURFACE AREA HERE %%%%%%%%%%
function [A] = surfaceArea(X,T)
    A = 0;
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE %%%%%%%%%
function [L] = cotLaplacian(X,T)
    nv = size(X,1);
    L = sparse(1:nv,1:nv,0);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE DIVIDED DIFFERENCES APPROXIMATION HERE %%%%%
function [G] = gradApprox(X,T)
    nv = size(X,1);
    G = zeros(nv,3);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%%

% ADD CODE TO COMPUTE THE BARYCENTRIC AREA VECTOR HERE %%%%%%%%%
function [M] = barycentricArea(X,T)
    nv = size(X,1);
    M = zeros(nv,1);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%

% ADD CODE TO COMPUTE POINTWISE MEAN CURVATURE HERE %%%%%%%%%%%
function [H] = meanCurvature(X,T)
    nv = size(X,1);
    H = zeros(nv,1);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%