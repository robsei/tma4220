clear all; clc; close all;

addpath('../grids');
addpath('../img');
addpath('../quadrature');

% Load image.
U = mat2gray(imread('cameraman.tif'));
U = imresize(U, 1);

% Parameters for anisotropic diffusion.
sigma = 0.00001;
tau = 0.005;
g = @(x) exp(-(x./2).^2);

% Order of quadrature.
Nq = 1;

% Grid generation.
[nx1, nx2] = size(U);
% [p, tri, edge] =  getSquareTri(nx1,nx2);
% u = U(:);
% [p, tri, edge] =  getSquareTri_swp(nx1, nx2, u);
[p, tri, edge, p_idx] =  getVarigrid(U, 0.001, 50);

% Image as vector.
U = U(sub2ind(size(U), p_idx(:,1), p_idx(:,2)));
u = U(:);

%dx1 = 1/(nx1-1);
%dx2 = 1/(nx2-1);

% Verify grid.
figure;
% imtrisurf(tri, p, u);
figure;
triplot(tri, p(:,1), p(:,2));

% Number of nodes.
n = size(p,1);
% Number of elements.
N = size(tri, 1);

% Initialize linear system.
M = sparse(n,n); % Mass matrix.
A_iso = sparse(n,n); % Isotropic stiffness matrix.
A_ani = sparse(n,n); % Anisotropic stiffness matrix.
b = zeros(n,1); % The right hand side.

disp('Started isotropic diffusion.');
tic;

% Each triangle is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tri(k,:);
    
    % Get the coordinated of the points.
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    
    % Compute phi's that land in this triangle.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    for i = 1:3
        for j = 1:3
            f_iso = @(x) dot(phi(1:2,i), phi(1:2,j));
            A_iso(nid(i),nid(j)) = A_iso(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,f_iso);
            
            f_mass = @(x) (phi(:,i)' * [x; 1]) * (phi(:,j)' * [x; 1]);
            M(nid(i),nid(j)) = M(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,f_mass);
            
            % g_rhs = @(x) U(nid(j)) * (phi(:,i)' * [x; 1]) * (phi(:,j)' * [x; 1]);
            % b(nid(i)) = b(nid(i)) +  quadrature2D(p1,p2,p3,Nq,g_rhs);
        end
    end
end
disp(['The loop took ' num2str(toc) 's.']);

u_iso = (M + sigma * A_iso) \ (M*u);

figure;
imtrisurf(tri, p, u_iso);

return

% Compute the diffusion coefficients.
d = zeros(N,2);

disp('Started computation of diffusion coefficients.');
tic;

% Each triangle is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tri(k,:);
    
    % Get the coordinated of the points.
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    
    % Compute phi's that land in this triangle.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    for i = 1:3
        d(k,:) = d(k,:) + u_iso(nid(i)) * phi(1:2,i)';
    end
end
disp(['The loop took ' num2str(toc) 's.']);

% Apply transfer function.
g_ani = g(sqrt(d(:,1).^2+d(:,2).^2));

% Visualize the diffusion coeffs
figure;
% plotdiffcoeff(nx1, nx2, g_ani);

disp('Started anisotropic diffusion.');
tic;

% Each triangle is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tri(k,:);
    
    % Get the coordinated of the points.
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    
    % Compute phi's that land in this triangle.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    for i = 1:3
        for j = 1:3
            f_ani = @(x) g_ani(k) * dot(phi(1:2,i), phi(1:2,j));
            A_ani(nid(i),nid(j)) = A_ani(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,f_ani);
        end
    end
end
disp(['The loop took ' num2str(toc) 's.']);
getVarigrid(U, 0.000001, 5000);
u_ani = (M + tau * A_ani) \ (M*u);

% Visualize the results
figure;
imtrisurf(tri, p, u_iso);

figure;
imtrisurf(tri, p, u_ani);