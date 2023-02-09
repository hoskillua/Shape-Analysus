% process triangle mesh
function mesh = MeshData(X,T)

if nargin == 0
    [X,T]=readOff('meshes/teddy.off'); 
end

mesh = [];
mesh.vertices = X;
mesh.triangles = double(T);

mesh.nv = size(X,1);
mesh.nf = size(T,1);

% compute triangle face normals
normalf = cross( mesh.vertices(mesh.triangles(:,2),:)'-mesh.vertices(mesh.triangles(:,1),:)', ...
                 mesh.vertices(mesh.triangles(:,3),:)'-mesh.vertices(mesh.triangles(:,1),:)' );
d = sqrt( sum(normalf.^2,1) ); d(d<eps)=1;
mesh.faceNormals = (normalf ./ repmat( d, 3,1 ))';

% compute triangle barycenters
mesh.faceCenters = (mesh.vertices(mesh.triangles(:,1),:)+mesh.vertices(mesh.triangles(:,2),:)+mesh.vertices(mesh.triangles(:,3),:))/3;

% compute mesh edges.
rawedges = mesh.triangles(:,[1 2 2 3 3 1])';
edges = reshape(rawedges,2,[])';
mesh.edges = unique(sort(edges,2),'rows');
mesh.ne = size(mesh.edges,1);

% compute triangle areas
x1 = mesh.vertices(mesh.triangles(:,1),1);
y1 = mesh.vertices(mesh.triangles(:,1),2);
z1 = mesh.vertices(mesh.triangles(:,1),3);
x2 = mesh.vertices(mesh.triangles(:,2),1);
y2 = mesh.vertices(mesh.triangles(:,2),2);
z2 = mesh.vertices(mesh.triangles(:,2),3);
x3 = mesh.vertices(mesh.triangles(:,3),1);
y3 = mesh.vertices(mesh.triangles(:,3),2);
z3 = mesh.vertices(mesh.triangles(:,3),3);
A = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
B = sqrt((x2-x3).^2+(y2-y3).^2+(z2-z3).^2);
C = sqrt((x3-x1).^2+(y3-y1).^2+(z3-z1).^2);
S = (A+B+C)/2;
mesh.triangleAreas = sqrt(S.*(S-A).*(S-B).*(S-C));

% compute edge lengths
mesh.edgeLengths = sqrt(sum((mesh.vertices(mesh.edges(:,1),:) - mesh.vertices(mesh.edges(:,2),:)).^2,2));

%{ 
Compute primal incidence matrix. This is an |E| x |V| matrix.
Row e has a +1, -1 pair at the vertex indices defining edge e.
The rest are zeros.
%}
mesh.primalIncidence = sparse(repmat(1:mesh.ne,1,2), mesh.edges(:), [ones(mesh.ne,1); -ones(mesh.ne,1)], mesh.ne, mesh.nv);

% compute area weighted vertex normals
T2V = sparse(repmat(1:mesh.nf,3,1),mesh.triangles',ones(numel(mesh.triangles),1),mesh.nf,mesh.nv);
mesh.vertNormals = T2V'*(mesh.faceNormals.*mesh.triangleAreas);
mesh.vertNormals = mesh.vertNormals ./ vecnorm(mesh.vertNormals,2,2);

% construct Cotan weights per edge
V2E = sparse(repmat(1:mesh.ne,2,1),mesh.edges',ones(numel(mesh.edges),1),mesh.ne,mesh.nv)';
T2E = T2V*V2E; T2E = T2E == 2;
[ii1, jj1] = find(T2E');
[ii2, jj2] = find(T2E(:,1:mesh.ne));
edges2triangles = sort(reshape(ii2,2,[])',2); 
v123 = mesh.triangles(edges2triangles(:,1),:);
v234 = mesh.triangles(edges2triangles(:,2),:);
v2 = mesh.edges(:,1);
v3 = mesh.edges(:,2);
[ii,jj] = find((v123 - v2).*(v123 - v3)); [~,perm] = sort(ii);
v1 = v123(sub2ind(size(v123), [1:mesh.ne]', jj(perm)));
[ii,jj] = find((v234 - v2).*(v234 - v3)); [~,perm] = sort(ii);
v4 = v234(sub2ind(size(v234), [1:mesh.ne]', jj(perm)));
edge2butterflyStencil = [v1 v2 v3 v4];
v1 = mesh.vertices(edge2butterflyStencil(:,1),:);
v2 = mesh.vertices(edge2butterflyStencil(:,2),:);
v3 = mesh.vertices(edge2butterflyStencil(:,3),:);
v4 = mesh.vertices(edge2butterflyStencil(:,4),:);
a1 = angleBetween(v2-v1,v3-v1);
a2 = angleBetween(v2-v4,v3-v4);
mesh.cotweights = cot(a1) + cot(a2);

% Construct cotan laplacian
mesh.cotLaplacian = mesh.primalIncidence' * diag(sparse(mesh.cotweights)) * mesh.primalIncidence /2;

% Construct barycentric vertex area. (sum of adjacent triangle areas divided by 3)
mesh.vertexWeights = accumarray(mesh.triangles(:),repmat(mesh.triangleAreas,3,1)/3,[mesh.nv,1]);
mesh.massMatrix = spdiags(mesh.vertexWeights, 0, mesh.nv, mesh.nv);

mesh.FtoV = mesh.massMatrix \ sparse(mesh.triangles, repmat((1:mesh.nf)', 1, 3), repmat(mesh.triangleAreas / 3, 1, 3), mesh.nv, mesh.nf);
mesh.VtoF = sparse(repmat((1:mesh.nf)', 1, 3), mesh.triangles, 1/3, mesh.nf, mesh.nv);

end

% computes vectorized angle between vector u and vector v
% u and v must both be Nx3 lists of vectors.
function z = angleBetween(u,v)
    u = u./vecnorm(u,2,2);
    v = v./vecnorm(v,2,2);
    z = acos(dot(u,v,2));
end
