function geodesicDistance()

addpath('utils/');

% load mesh
meshfile = 'meshes/moomoo.off';
[X,T] = readOff(meshfile);
data = MeshData(X,T);

% build signal: delta heat source at highest Y axis point. This defines the source vertex.
dt = sqrt(mean(data.edgeLengths));
signal = zeros(data.nv, 1);
[~, maxind] = max(data.vertices(:, 2));
signal(maxind) = 1;

%% PROBLEM 2 - THE HEAT METHOD
%%% YOUR CODE HERE - short-time heat diffusion
phi = zeros(data.nv, 1);

[Div, Grad] = getDivGrad(data);

%%% YOUR CODE HERE - compute the normalized gradient of phi
qvec = zeros(data.nf, 3);

%%% YOUR CODE HERE - compute divergence of normalized gradient
divVals = zeros(data.nv, 1);

%%% YOUR CODE HERE - final geodesic distances
dist = zeros(data.nv, 1);
%%% END HOMEWORK PROBLEM

%% Visualization
figure;
tlt = tiledlayout(1, 3);
tlt.TileSpacing = 'compact';
tlt.Padding = 'compact';

nexttile(1);
trisurf(T, X(:, 1), X(:, 2), X(:, 3), phi, 'LineStyle', 'none'); view(3); axis image vis3d off; shading interp;
hold on;
quiver3(data.faceCenters(:,1), data.faceCenters(:,2), data.faceCenters(:,3), ...
        qvec(:,1), qvec(:,2), qvec(:,3), 'r');
view(3); axis image vis3d off; shading interp;
title('Heat kernel & normalized gradient');

nexttile(2);
trisurf(T, X(:, 1), X(:, 2), X(:, 3), divVals, 'LineStyle', 'none'); view(3); axis image vis3d off; shading interp;
title('Divergences');

ax = nexttile(3);
trisurf(T, X(:, 1), X(:, 2), X(:, 3), dist, 'LineStyle', 'none'); view(3); axis image vis3d off; shading interp; colormap(ax, parula(20));
title('Geodesic distances');

end
