function [p, tri, edge] = getGradgrid(U, sigma, g)

% Blur image to remove noise.
U = imgaussfilt(U, sigma);

% Compute vertical and horizontal derivatives.
Diff_x1 = conv2(U, [1;-1;0], 'same');
Diff_x2 = conv2(U, [1,-1,0], 'same');

Diff_x1(end,:) = 0;
Diff_x2(:,end) = 0;

% Get normalized norm of gradient.
Diff = sqrt(Diff_x1.^2 + Diff_x2.^2);

d0 = min(Diff(:));
d1 = max(Diff(:));

Diff = (Diff - d0) / (d1 - d0);

% Apply transfer function.
Diff = g(Diff);

% Compute mask.
% Points that will be used as nodes are set to one.
mask = (rand(size(U)) < Diff);
mask(1,1) = 1;
mask(1,end) = 1;
mask(end,1) = 1;
mask(end,end) = 1;

[p, tri, edge] = maskDelaunay(mask, size(U));

end

