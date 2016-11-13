clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(1, 0.05);

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 0.8;
g = @(x) exp(-(12*x).^2);

% Grid generation.
[p, tri, edge] = getGradgrid(U, 2, @(x) x);
u = maskApply(U, p);

% Run PM.
[u_ani, ~, ~, ~, ~, ~] = pm_tri_pre(u, tri, p, edge, sigma, tau, g);

% Plot solution.
figure;
subplot(2,2,1);
imshow(U);
subplot(2,2,2);
imtriplot(tri, p, edge);
subplot(2,2,3);
imtrisurf(tri, p, u_ani);
subplot(2,2,4);
iminterpsurf(tri, p, u_ani, size(U));

% Print PSNR.
disp(['PSNR is ' num2str(psnr(u_ani, maskApply(U_true, p)))]);

% Write results to file.
print('../out/report_gradgrid1','-dpng');

% Parameters for anisotropic diffusion.
sigma = 0.5;
tau = 1.4;
g = @(x) exp(-(12*x).^2);

% Run PM.
[u_ani, ~, ~, ~, ~, ~] = pm_tri_pre(u, tri, p, edge, sigma, tau, g);

% Plot solution.
figure;
subplot(2,2,1);
imshow(U);
subplot(2,2,2);
imtriplot(tri, p, edge);
subplot(2,2,3);
imtrisurf(tri, p, u_ani);
subplot(2,2,4);
iminterpsurf(tri, p, u_ani, size(U));

% Print PSNR.
disp(['PSNR is ' num2str(psnr(u_ani, maskApply(U_true, p)))]);

% Write results to file.
print('../out/report_gradgrid2','-dpng');