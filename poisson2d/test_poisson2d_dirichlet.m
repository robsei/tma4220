clear all; clc; close all;

% Load needed functions.
addpath('../grids');

% Boundary conditions.
f = @(X) 16 .* pi^2 .* X(1,:) .* X(2,:) .* (X(1,:).^2 + X(2,:).^2) .* sin(2.*pi.*(X(1,:).^2 + X(2,:).^2)) - 24.*X(1,:).*X(2,:).*pi.*cos(2.*pi.*(X(1,:).^2 + X(2,:).^2));

% Grid parameters.
ns = [50, 500, 1000];
n_max = length(ns);
theta = 3/2 * pi;

% Solve problem for different n.
for i = 1:n_max
    % Construct the grid.
    [p, tri, edge] = getSlice(ns(i), theta);
    
    % Solve problem
    u_fem = poisson2d_dirichlet(p, tri, edge, f);
    
    % Plot FEM solution.
    subplot(n_max, 2, i*2-1);
    trimesh(tri, p(:,1), p(:,2), u_fem);

    % Plot analytic solution.
    subplot(n_max, 2, i*2);
    u_ana = p(:,1) .* p(:,2) .* sin(2.*pi.*(p(:,1).^2 + p(:,2).^2));
    trimesh(tri, p(:,1), p(:,2), u_ana);
end

% Write results to file.
print('../out/test_poisson2d_dirichlet','-dpng');