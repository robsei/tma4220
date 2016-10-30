function [p, tri, edge] = maskDelaunay(mask, nx)

x1 = 1:nx(1);
x2 = 1:nx(2);
[X2, X1] = meshgrid(x1, x2);
x(:,1) =  X1(mask);
x(:,2) =  X2(mask);

dt = delaunayTriangulation(x);
p = dt.Points;
tri = dt.ConnectivityList;
edge = freeBoundary(dt);

end
