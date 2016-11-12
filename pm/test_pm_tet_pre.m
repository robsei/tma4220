clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
U = loadCameraman(0.3, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 1;
g = @(x) exp(-(8*x).^2);
%g = @(x) exp(-(5*x).^2);
%g = @(x) -tanh(5*x)+1;
%g = @(x) ones(size(x));
%g = @(x) 1-x;

% Grid generation.
[nx1, nx2] = size(U);
[p, tet, edge] = getSquareTet(nx1, nx2);

% Run PM.
[U_iso, ~, U_ani] = pm_tet_pre(U, sigma, tau, g);

% Plot solution.
subplot(1,3,1);
imshow(U);
subplot(1,3,2);
imshow(U_iso);
subplot(1,3,3);
imshow(U_ani);

print('../out/test_pm_tet_pre', '-dpng');
