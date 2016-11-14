clear all; clc; close all;

addpath('../grids');
addpath('../img');

% Load image.
[U, U_true] = loadCameraman(1, 0.1);

% Parameters for anisotropic diffusion.
sigma = 0.5;
taus = exp(-9:0.05:9);
imax = length(taus);
g = @(x) exp(-(12*x).^2);

% Grid generation.
[p, tri, edge] = getVarigrid(U, 2, 0.0001, 0.2, 4, 50);
u = maskApply(U, p);

% Allocate psnrs.
psnrs = zeros(imax,1);

% Prepare video output.
v1 = VideoWriter('../out/video_taus1.avi');
v2 = VideoWriter('../out/video_taus2.avi');
v3 = VideoWriter('../out/video_taus3.avi');
open(v1);
open(v2);
open(v3);

% Run PM.
[~, A_iso, M, ~, ~, A_ani] = pm_tri_pre(u, tri, p, edge, sigma, 0, g);

% Precompute rhs.
b = M*u;

for i=1:imax
    % Solve linear systems.
    u_iso = (M + taus(i) * A_iso) \ b;
    u_ani = (M + taus(i) * A_ani) \ b;

    % Plot solutions.
    figure;
    iminterpsurf(tri, p, u_iso, size(U));

    % Write results to video.
    writeVideo(v1, getframe());
    
    % Plot solutions.
    figure;
    iminterpsurf(tri, p, u_ani, size(U));

    % Write results to video.
    writeVideo(v2, getframe());

    % Compute psnr.
    psnrs(i) = impsnr(p, u_ani, U_true);

    % Plot solutions.
    figure;
    subplot(1,2,1);
    iminterpsurf(tri, p, u_ani, size(U));
    subplot(1,2,2);
    semilogx(taus, psnrs);
    axis square;

    % Write results to video.
    writeVideo(v3, getframe(gcf()));
    
    disp(['Frame ' num2str(i) ' of ' num2str(imax)]);
    close all;
end

close(v1);
close(v2);
close(v3);
