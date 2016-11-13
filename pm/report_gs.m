clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(1, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 4;
K = [-100 -10 -1 -0.1];

for i=1:length(K)
    % Grid generation.
    [p, tri, edge] = getVarigrid(U, 2, 0.0001, 0.2, 4, 50);
    u = maskApply(U, p);

    % Run PM.
    g = @(x) exp(-(K(i)*x).^2);
    [u_ani, ~, ~, ~, ~, ~] = pm_tri_pre(u, tri, p, edge, sigma, tau, g);
    
    % Plot solution.
    subplot(2,length(K)/2,i);
    iminterpsurf(tri, p, u_ani, size(U));

    % Print PSNR.
    disp(['PSNR is ' num2str(impsnr(p, u_ani, U_true))]);
end

% Write results to file.
print('../out/report_gs','-dpng');