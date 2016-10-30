clear all; clc; close all;

% Load needed functions.
addpath('../quadrature');
addpath('../grids');

% Boundary conditions.
f = @(X) 16 .* pi^2 .* X(1,:) .* X(2,:) .* (X(1,:).^2 + X(2,:).^2) .* sin(2.*pi.*(X(1,:).^2 + X(2,:).^2)) - 24.*X(1,:).*X(2,:).*pi.*cos(2.*pi.*(X(1,:).^2 + X(2,:).^2));
g_1active = @(X) -X * sin(2*pi*X.^2);
g_2active = @(X) X * sin(2*pi*X.^2);

% Grid construction.
n = 500;
theta = 3/2 * pi;
[p, tri, edge] = getSlice(n, theta);

N = size(tri, 1);

A = zeros(n,n);
b = zeros(n,1);

triplot(tri, p(:,1), p(:,2));

% Each triangle is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tri(k,:);
    
    % Compute phi's that land in this triangle.
    % Create lhs of the linear system.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    % Is the triangle on the Neumann boundary?
    is_on_edge = (abs(p(nid,1)) < eps & p(nid,2) < eps) | (p(nid,1) > eps & abs(p(nid,2)) < eps);
    is_nm_edge = sum(is_on_edge) == 2;
    if is_nm_edge
        pt_nm_edge = p(nid(is_on_edge),:);
        active_dim = sum(abs(pt_nm_edge),1) > (2*eps);
    end
    
    
    for i = 1:3
        p1 = p(nid(1),:);
        p2 = p(nid(2),:);
        p3 = p(nid(3),:);
        Nq = 1;
        for j = 1:3N = size(tri, 1);

A = zeros(n,n);
b = zeros(n,1);

triplot(tri, p(:,1), p(:,2));

% Each triangle is visited once.
for k = 1:N
    % Get node ids to identify points and basis functions.
    nid = tri(k,:);
    
    % Compute phi's that land in this triangle.
    % Create lhs of the linear system.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    % Is the triangle on the Neumann boundary?
    is_on_edge = (abs(p(nid,1)) < eps & p(nid,2) < eps) | (p(nid,1) > eps & abs(p(nid,2)) < eps);
    is_nm_edge = sum(is_on_edge) == 2;
    if is_nm_edge
        pt_nm_edge = p(nid(is_on_edge),:);
        active_dim = sum(abs(pt_nm_edge),1) > (2*eps);
    end
    
    
            g1 = @(x) dot(phi(1:2,i), phi(1:2,j)) * ones(1, size(x,2));
            A(nid(i),nid(j)) = A(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,g1);
        end
        g2 = @(x) f(x) * phi(:,i)' * [x; ones(1, size(x,2))];
        g3_1active = @(x) g_1active(x) * phi(:,i)' * [x; 0; 1];
        g3_2active = @(x) g_2active(x) * phi(:,i)' * [0; x; 1];
        b(nid(i)) = b(nid(i)) + quadrature2D(p1,p2,p3,Nq,g2);
        if is_nm_edge
            if active_dim(1)
                b(nid(i)) = b(nid(i)) + quadrature1D(pt_nm_edge(1,active_dim),pt_nm_edge(2,active_dim),Nq,g3_1active);
            elseif active_dim(2)
                b(nid(i)) = b(nid(i)) + quadrature1D(pt_nm_edge(1,active_dim),pt_nm_edge(2,active_dim),Nq,g3_2active);
            end
        end
    end
end

A(edge(:,1),:) = [];
A(:,edge(:,1)) = [];
b(edge(:,1)) = [];
u = A \ b;

U = zeros(n,1);
idx = 1:n;
idx(edge(:,1)) = [];
U(idx) = u;

subplot(1,2,1);
trimesh(tri, p(:,1), p(:,2), U);

% Verify the result.
u_true = p(:,1) .* p(:,2) .* sin(2.*pi.*(p(:,1).^2 + p(:,2).^2));
subplot(1,2,2);
trimesh(tri, p(:,1), p(:,2), u_true);
disp(['l2-error of solution is ', num2str(norm(U - u_true))]);