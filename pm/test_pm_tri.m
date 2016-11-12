clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
U = loadCameraman(0.2, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 1.5;
g = @(x) exp(-(10*x).^2);
%g = @(x) exp(-(5*x).^2);
%g = @(x) -tanh(5*x)+1;
%g = @(x) ones(size(x));
%g = @(x) 1-x;

% Grid generation.
% [p, tri, edge] = getGradgrid(U, 2, @(x) power(x,1));
% [p, tri, edge] = getVarigrid(U, 2, 0.0001, 0.2, 5, 50);
[p, tri, edge] = getSquareTri_swp(U); 
u = maskApply(U, p);

% Run PM.
[u_ani, ~, ~, u_iso, ~, ~] = pm_tri(u, tri, p, edge, sigma, tau, g);

% Plot solution.
subplot(2,2,1);
imshow(U);
subplot(2,2,2);
imtriplot(tri, p, edge);
subplot(2,2,3);
iminterpsurf(tri, p, u_iso, size(U));
subplot(2,2,4);
iminterpsurf(tri, p, u_ani, size(U));

% Write results to file.
print('../out/test_pm_tri','-dpng');
