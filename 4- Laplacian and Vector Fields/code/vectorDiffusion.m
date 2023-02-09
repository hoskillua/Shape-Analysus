function vectorDiffusion()

addpath('utils/');

% load mesh
meshfile = 'meshes/moomoo.off';
[X,T] = readOff(meshfile);
data = MeshData(X,T);
CL = buildConnectionLaplacian(data);

% build vector valued signal: delta vector heat source at highest Y axis point.
dt = sqrt(mean(data.edgeLengths));
signal = zeros(3*data.nv,1);
[~,maxind] = max(data.vertices(:,2));
t1 = cross(data.vertNormals(maxind,:),randn(3,1)'); t1 = t1/norm(t1);
signal((maxind-1)*3 + (1:3)) = t1;

%% PROBLEM 3(b) - VECTOR HEAT METHOD
% todo: compute diffused signal for various dt and normalize
phi = zeros(data.nv, 3);
%%% END HOMEWORK PROBLEM

%% Visualization
figure; trisurf(T, X(:, 1), X(:, 2), X(:, 3), 'FaceColor', 'k', 'LineStyle', 'none'); hold on;
q2 = quiver3(X(:,1),X(:,2),X(:,3),phi(:,1),phi(:,2),phi(:,3), 'r', 'LineWidth', 1);
q1 = quiver3(X(maxind,1),X(maxind,2),X(maxind,3),t1(1),t1(2),t1(3), 'c', 'LineWidth', 2, 'AutoScaleFactor', 3);
view(3); axis image vis3d off;
title('parallel transported vector field');
legend([q1, q2], 'Source', 'Parallel transport');

end

%% PROBLEM 3(a)
%%% YOUR CODE HERE - implement the connection laplacian operator
% data: mesh data structure
% CL: |3V| x |3V| "sparse connection laplacian operator"
function CL = buildConnectionLaplacian(data)
    CL = sparse(3 * data.nv, 3 * data.nv);
end
%%% END HOMEWORK PROBLEM

