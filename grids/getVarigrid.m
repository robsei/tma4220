function [p, tri, edge] = getVarigrid(U, sigma, minvar, cvc, alldepth, maxdepth)

% Blur image to remove noise.
if sigma > 0
    U = imgaussfilt(U, sigma);
end

[nx1, nx2] = size(U);
idx1 = 1:nx1;
idx2 = 1:nx2;
depth = 0;

p = getVarigrid_recursion(U, idx1, idx2, minvar, cvc, depth, alldepth, maxdepth);

p(end+1,:) = [1 1];
p(end+1,:) = [1 nx2];
p(end+1,:) = [nx1 1];
p(end+1,:) = [nx1 nx2];

dt = delaunayTriangulation(p);
p = dt.Points;
tri = dt.ConnectivityList;
edge = freeBoundary(dt);

end

