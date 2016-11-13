clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(1, 0.05);

% Parameters for grid generation.
powers = 1:3;

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 1.5;
g = @(x) exp(-(12*x).^2);

for i=1:length(powers)
    % Grid generation.
    transfer = @(X) power(X,powers(i));
    [p, tri, edge] = getGradgrid(U, 2, transfer);
    u = maskApply(U, p);

    % Run PM.
    [u_ani, ~, ~, ~, ~, ~] = pm_tri_pre(u, tri, p, edge, sigma, tau, g);
    
    % Plot solution.
    subplot(2,length(powers),i);
    imtriplot(tri, p, edge);
    subplot(2,length(powers),i+length(powers));
    iminterpsurf(tri, p, u_ani, size(U));

    % Print PSNR.
    disp(['PSNR is ' num2str(impsnr(p, u_ani, U_true))]);
end

% Write results to file.
print('../out/report_transfers','-dpng');