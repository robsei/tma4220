clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(0.5, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 0.8;
g = @(x) exp(-(12*x).^2);

% Grid generation.
[p, tri, edge] = getSquareTri(size(U)); 
u = maskApply(U, p);

% Run PM.
[u_ani, ~, ~, ~, ~, ~] = pm_tri(u, tri, p, edge, sigma, tau, g);

% Plot solution.
subplot(2,2,1);
imshow(U);
subplot(2,2,2);
imtriplot(tri, p, edge);
subplot(2,2,3);
imtrisurf(tri, p, u_ani);
subplot(2,2,4);
iminterpsurf(tri, p, u_ani, size(U));

% Print PSNR.
disp(['PSNR is ' num2str(impsnr(p, u_ani, U_true))]);

% Write results to file.
print('../out/report_firsttry','-dpng');