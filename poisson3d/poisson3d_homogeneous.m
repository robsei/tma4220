clear all; clc; close all;

% Load needed functions.
addpath('../quadrature');
addpath('../grids');

% PDE right hand side.
f = @(X) 12 .* pi^2 .* sin(2 .* pi .* X(1)) .* sin(2 .* pi .* X(2)) .* sin(2 .* pi .* X(3));

% Grid construction.
n = 15;
[p, tet, edge] = getBox(n);

N = size(tet, 1);

A = zeros(n^3,n^3);
b = zeros(n^3,1);

% Each thetrahedron is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tet(k,:);
    
    % Compute phi's that land in this triangle.
    % Create lhs of the linear system.
    lhs = [p(nid,:), ones(4,1)];
    phi = zeros(4,4);
    phi(:,1) = lhs \ [1; 0; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0; 0];
    phi(:,3) = lhs \ [0; 0; 1; 0];
    phi(:,4) = lhs \ [0; 0; 0; 1];
    
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    p4 = p(nid(4),:);
    for i = 1:4
        Nq = 1;
        for j = 1:4
            g1 = @(x) dot(phi(1:3,i), phi(1:3,j));
            A(nid(i),nid(j)) = A(nid(i),nid(j)) + quadrature3D(p1,p2,p3,p4,Nq,g1);
        end
        g2 = @(x) f(x) * phi(:,i)' * [x; 1];
        b(nid(i)) = b(nid(i)) +  quadrature3D(p1,p2,p3,p4,Nq,g2);
    end
end

all_edge_nodes = unique(edge);
A(all_edge_nodes,:) = [];
A(:,all_edge_nodes) = [];
b(all_edge_nodes) = [];
u = A \ b;

U = zeros(n^3,1);
idx = 1:(n^3);
idx(all_edge_nodes) = [];
U(idx) = u;

% l2 error
u_true = sin(2 .* pi .* p(:,1)) .* sin(2 .* pi .* p(:,2)) .* sin(2 .* pi .* p(:,3));
disp(['l2-error of solution is ', num2str(norm(U - u_true))]);

% New verification procedure.
z_levels = unique(p(:,3));
for z_level = z_levels'
    z_idx = p(:,3) == z_level;
    X = p(z_idx,1);
    Y = p(z_idx,2);
    Z = p(z_idx,3);
    V_comp = U(z_idx);
    V_true = sin(2 .* pi .* X) .* sin(2 .* pi .* Y) .* sin(2 .* pi .* Z);
    
    subplot(1,2,1);
    plot3(X,Y,V_comp,'.');
    subplot(1,2,2);
    plot3(X,Y,V_true,'.');
    pause;
end
% 
% % Old verification precedure.
% u_interp = scatteredInterpolant(p,U);
% lift = 0.0;
% [X,Y] = meshgrid(0:0.1:1);
% v = u_interp(X(:), Y(:), lift * ones(size(X(:),1),1));
% subplot(1,2,1);
% plot3(X(:), Y(:), v, '.');
% 
% % Verify the result.
% u_true = @(X,Y,Z) sin(2 .* pi .* X) .* sin(2 .* pi .* Y) .* sin(2 .* pi .* Z);
% v = u_true(X(:), Y(:), lift * ones(size(X(:),1),1));
% subplot(1,2,2);
% plot3(X(:), Y(:), v, '.');