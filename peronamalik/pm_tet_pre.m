clear all; clc; close all;

addpath('../grids');

% Load and prepare test data.
U = mat2gray(imread('cameraman.tif'));
% U = U(:,:,1);
% U = U + 0.1 * randn(size(U));
U = imresize(U,0.4);

% Grid generation.
[nx1, nx2] = size(U);
dx1 = 1/(nx1-1);
dx2 = 1/(nx2-1);
[p, tet, edge] = getSquareTet(nx1, nx2);

% Settings for PM equation.
sigma = 0.0001;
tau = 0.0025;
g = @(x) max(exp(-(x./2).^2), 0);

% Number of nodes.
n = nx1 * nx2;
% Number of elements.
N = size(tet, 1);

% Allocate memory.
M = sparse(n,n); % Mass matrix.
A_iso = sparse(n,n); % Isotropic stiffness matrix.
A_ani = sparse(n,n); % Anisotropic stiffness matrix.
b = zeros(n,1); % The right hand side.

% Loop for isotropic diffusion.
% Each triangle is visited once.
for i = 0:nx1-2
    for j = 0:nx2-2
        A_iso(i + j*nx1 + 1,i + j*nx1 + 1) = A_iso(i + j*nx1 + 1,i + j*nx1 + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_iso(i + j*nx1 + 1,i + nx1*(j + 1) + 1) = A_iso(i + j*nx1 + 1,i + nx1*(j + 1) + 1) + -dx1*dx2*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_iso(i + j*nx1 + 1,i + nx1*(j + 1) + 2) = A_iso(i + j*nx1 + 1,i + nx1*(j + 1) + 2) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + j*nx1 + 1,i + j*nx1 + 2) = A_iso(i + j*nx1 + 1,i + j*nx1 + 2) + -dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_iso(i + nx1*(j + 1) + 1,i + j*nx1 + 1) = A_iso(i + nx1*(j + 1) + 1,i + j*nx1 + 1) + -dx1*dx2*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_iso(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) = A_iso(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) + dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_iso(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) = A_iso(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) + -dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A_iso(i + nx1*(j + 1) + 1,i + j*nx1 + 2) = A_iso(i + nx1*(j + 1) + 1,i + j*nx1 + 2) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + nx1*(j + 1) + 2,i + j*nx1 + 1) = A_iso(i + nx1*(j + 1) + 2,i + j*nx1 + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) = A_iso(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) + -dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A_iso(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) = A_iso(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) + dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + (j/dx1 - (j + 1/2)/dx1)^2);
        A_iso(i + nx1*(j + 1) + 2,i + j*nx1 + 2) = A_iso(i + nx1*(j + 1) + 2,i + j*nx1 + 2) + -dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + j*nx1 + 2,i + j*nx1 + 1) = A_iso(i + j*nx1 + 2,i + j*nx1 + 1) + -dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_iso(i + j*nx1 + 2,i + nx1*(j + 1) + 1) = A_iso(i + j*nx1 + 2,i + nx1*(j + 1) + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + j*nx1 + 2,i + nx1*(j + 1) + 2) = A_iso(i + j*nx1 + 2,i + nx1*(j + 1) + 2) + -dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_iso(i + j*nx1 + 2,i + j*nx1 + 2) = A_iso(i + j*nx1 + 2,i + j*nx1 + 2) + dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        M(i + j*nx1 + 1,i + j*nx1 + 1) = M(i + j*nx1 + 1,i + j*nx1 + 1) + dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
        M(i + j*nx1 + 1,i + nx1*(j + 1) + 1) = M(i + j*nx1 + 1,i + nx1*(j + 1) + 1) + -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + nx1*(j + 1) + 2) = M(i + j*nx1 + 1,i + nx1*(j + 1) + 2) + dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + j*nx1 + 2) = M(i + j*nx1 + 1,i + j*nx1 + 2) + -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 1,i + j*nx1 + 1) = M(i + nx1*(j + 1) + 1,i + j*nx1 + 1) + -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) = M(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) + dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))^2;
        M(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) = M(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) + -dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 1,i + j*nx1 + 2) = M(i + nx1*(j + 1) + 1,i + j*nx1 + 2) + dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 2,i + j*nx1 + 1) = M(i + nx1*(j + 1) + 2,i + j*nx1 + 1) + dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) = M(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) + -dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) = M(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) + dx1*dx2*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2))^2;
        M(i + nx1*(j + 1) + 2,i + j*nx1 + 2) = M(i + nx1*(j + 1) + 2,i + j*nx1 + 2) + -dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 2,i + j*nx1 + 1) = M(i + j*nx1 + 2,i + j*nx1 + 1) + -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 2,i + nx1*(j + 1) + 1) = M(i + j*nx1 + 2,i + nx1*(j + 1) + 1) + dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 2,i + nx1*(j + 1) + 2) = M(i + j*nx1 + 2,i + nx1*(j + 1) + 2) + -dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 2,i + j*nx1 + 2) = M(i + j*nx1 + 2,i + j*nx1 + 2) + dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
    end
end

% Try RHS like in lecture.
b = U(:);

% Treat boundary conditions.
A_iso(edge(:,1),:) = [];
A_iso(:,edge(:,1)) = [];
M(edge(:,1),:) = [];
M(:,edge(:,1)) = [];
b(edge(:,1)) = [];

% Solve linear system.
u_iso = (M + sigma * A_iso) \ (M*b);

% Fill in boundary condition.
u_iso_bd = zeros(n,1);
idx = 1:n;
idx(edge(:,1)) = [];
u_iso_bd(idx) = u_iso;

% u_iso_bd = u_iso;

figure;
subplot(1,2,1);
imshow(U);
subplot(1,2,2);
imshow(reshape(u_iso_bd,nx1,nx2));

d = zeros(n, 2);

for i = 0:nx1-2
    for j = 0:nx2-2
        d(i + j*nx1 + 1,1) = d(i + j*nx1 + 1,1) + -u_iso_bd(i + j*nx1 + 1)*((j + 1)/dx1 - (j + 1/2)/dx1);
        d(i + j*nx1 + 1,2) = d(i + j*nx1 + 1,2) + -u_iso_bd(i + j*nx1 + 1)*((i + 1)/dx2 - (i + 1/2)/dx2);
        d(i + j*nx1 + 1,1) = d(i + j*nx1 + 1,1) + u_iso_bd(i + nx1*(j + 1) + 1)*(j/dx1 - (j + 1/2)/dx1);
        d(i + j*nx1 + 1,2) = d(i + j*nx1 + 1,2) + u_iso_bd(i + nx1*(j + 1) + 1)*((i + 1)/dx2 - (i + 1/2)/dx2);
        d(i + j*nx1 + 1,1) = d(i + j*nx1 + 1,1) + -u_iso_bd(i + nx1*(j + 1) + 2)*(j/dx1 - (j + 1/2)/dx1);
        d(i + j*nx1 + 1,2) = d(i + j*nx1 + 1,2) + -u_iso_bd(i + nx1*(j + 1) + 2)*(i/dx2 - (i + 1/2)/dx2);
        d(i + j*nx1 + 1,1) = d(i + j*nx1 + 1,1) + u_iso_bd(i + j*nx1 + 2)*((j + 1)/dx1 - (j + 1/2)/dx1);
        d(i + j*nx1 + 1,2) = d(i + j*nx1 + 1,2) + u_iso_bd(i + j*nx1 + 2)*(i/dx2 - (i + 1/2)/dx2);
    end
end

figure;
subplot(1,3,1);
imagesc(reshape(d(:,1),nx1,nx2))
subplot(1,3,2);
imagesc(reshape(d(:,2),nx1,nx2))

g_ani = max(g(sqrt(d(:,1).^2+d(:,2).^2)), 0.01);


subplot(1,3,3);
imagesc(reshape(g_ani,nx1,nx2));

for i = 0:nx1-2
    for j = 0:nx2-2
        A_ani(i + j*nx1 + 1,i + j*nx1 + 1) = A_ani(i + j*nx1 + 1,i + j*nx1 + 1) + dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_ani(i + j*nx1 + 1,i + nx1*(j + 1) + 1) = A_ani(i + j*nx1 + 1,i + nx1*(j + 1) + 1) + -dx1*dx2*g_ani(i + j*nx1 + 1)*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_ani(i + j*nx1 + 1,i + nx1*(j + 1) + 2) = A_ani(i + j*nx1 + 1,i + nx1*(j + 1) + 2) + dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + j*nx1 + 1,i + j*nx1 + 2) = A_ani(i + j*nx1 + 1,i + j*nx1 + 2) + -dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_ani(i + nx1*(j + 1) + 1,i + j*nx1 + 1) = A_ani(i + nx1*(j + 1) + 1,i + j*nx1 + 1) + -dx1*dx2*g_ani(i + j*nx1 + 1)*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_ani(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) = A_ani(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) + dx1*dx2*g_ani(i + j*nx1 + 1)*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A_ani(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) = A_ani(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) + -dx1*dx2*g_ani(i + j*nx1 + 1)*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A_ani(i + nx1*(j + 1) + 1,i + j*nx1 + 2) = A_ani(i + nx1*(j + 1) + 1,i + j*nx1 + 2) + dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + nx1*(j + 1) + 2,i + j*nx1 + 1) = A_ani(i + nx1*(j + 1) + 2,i + j*nx1 + 1) + dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) = A_ani(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) + -dx1*dx2*g_ani(i + j*nx1 + 1)*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A_ani(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) = A_ani(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) + dx1*dx2*g_ani(i + j*nx1 + 1)*((i/dx2 - (i + 1/2)/dx2)^2 + (j/dx1 - (j + 1/2)/dx1)^2);
        A_ani(i + nx1*(j + 1) + 2,i + j*nx1 + 2) = A_ani(i + nx1*(j + 1) + 2,i + j*nx1 + 2) + -dx1*dx2*g_ani(i + j*nx1 + 1)*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + j*nx1 + 2,i + j*nx1 + 1) = A_ani(i + j*nx1 + 2,i + j*nx1 + 1) + -dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A_ani(i + j*nx1 + 2,i + nx1*(j + 1) + 1) = A_ani(i + j*nx1 + 2,i + nx1*(j + 1) + 1) + dx1*dx2*g_ani(i + j*nx1 + 1)*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + j*nx1 + 2,i + nx1*(j + 1) + 2) = A_ani(i + j*nx1 + 2,i + nx1*(j + 1) + 2) + -dx1*dx2*g_ani(i + j*nx1 + 1)*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A_ani(i + j*nx1 + 2,i + j*nx1 + 2) = A_ani(i + j*nx1 + 2,i + j*nx1 + 2) + dx1*dx2*g_ani(i + j*nx1 + 1)*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
    end
end

% Treat boundary conditions.
A_ani(edge(:,1),:) = [];
A_ani(:,edge(:,1)) = [];

u_ani = (M + tau * A_ani) \ (M*b);

% Fill in boundary condition.
u_ani_bd = zeros(n,1);
idx = 1:n;
idx(edge(:,1)) = [];
u_ani_bd(idx) = u_ani;

% u_ani_bd = u_ani;

figure;
subplot(1,2,1);
imshow(U);
subplot(1,2,2);
imshow(reshape(u_ani_bd,nx1,nx2));
