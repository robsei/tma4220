clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(1, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
taus = [0, 0.5, 1, 5, 10, 50, 100, 500, 1000, 5000];
g = @(x) exp(-(12*x).^2);

% Grid generation.
[p, tri, edge] = getVarigrid(U, 2, 0.0001, 0.2, 4, 50);
u = maskApply(U, p);

% Run PM.
[~, ~, M, ~, ~, A_ani] = pm_tri_pre(u, tri, p, edge, sigma, 0, g);

for i=1:length(taus)
    % Solve linear system.
    u_ani = (M + taus(i) * A_ani) \ (M*u);
    
    % Plot solution.
    figure;
    iminterpsurf(tri, p, u_ani, size(U));

    % Print PSNR.
    disp(['PSNR is ' num2str(impsnr(p, u_ani, U_true))]);

    % Write results to file.
    print(['../out/report_taus', num2str(i, '%.2d')],'-dpng');
end