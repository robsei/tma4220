function [X, tri, edge, P] = getGradgrid(U, sigma, g)

U = imgaussfilt(U, sigma);

[nx1, nx2] = size(U);
dx1 = 1/(nx1-1);
dx2 = 1/(nx2-1);

Diff_x1 = conv2(U, [1;-1;0], 'same');
Diff_x2 = conv2(U, [1,-1,0], 'same');

Diff_x1(end,:) = 0;
Diff_x2(:,end) = 0;

Diff = sqrt(Diff_x1.^2 + Diff_x2.^2);

d0 = min(Diff(:));
d1 = max(Diff(:));

Diff = (Diff - d0) / (d1 - d0);

Diff = g(Diff);

mask = (rand(size(U)) < Diff);
mask(1,1) = 1;
mask(1,end) = 1;
mask(end,1) = 1;
mask(end,end) = 1;

P1 = 1:nx1;
P2 = 1:nx2;
[P2, P1] = meshgrid(P1, P2);
P(:,1) =  P1(mask);
P(:,2) =  P2(mask);

x1 = linspace(0,1,nx1);
x2 = linspace(0,1,nx2);
[X2, X1] = meshgrid(x1, x2);

X(:,1) =  X1(mask);
X(:,2) =  X2(mask);

dt = delaunayTriangulation(X);
tri = dt.ConnectivityList;
edge = freeBoundary(dt);

end

