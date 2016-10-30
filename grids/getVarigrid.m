function [p, tri, edge, p_idx] = getVarigrid(U, minvar, maxdepth)

[nx1, nx2] = size(U);

idx1 = 1:nx1;
idx2 = 1:nx2;

depth = 0;

p_idx = getVarigrid_recursion(U, idx1, idx2, minvar, depth, maxdepth);

p_idx(end+1,:) = [1 1];
p_idx(end+1,:) = [1 nx2];
p_idx(end+1,:) = [nx1 1];
p_idx(end+1,:) = [nx1 nx2];

dx1 = 1/(nx1-1);
dx2 = 1/(nx2-1);

p = [(p_idx(:,1)-1)*dx1 (p_idx(:,2)-1)*dx2];

dt = delaunayTriangulation(p);
% p = dt.Points;
tri = dt.ConnectivityList;
edge = freeBoundary(dt);

% p_idx = [min(round(p(:,1)*nx1)+1,nx1) min(round(p(:,1)*nx1)+1,nx2)];

end

