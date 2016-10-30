function [u_iso, g_ani, u_ani] = pm_tri_pre(u, tri, p, edge, sigma, tau, g)

% Number of nodes.
n = size(p,1);
% Number of elements.
N = size(tri, 1);

% Precompute basis functions and quadrature for later use.
phi = zeros(3, 3, N);
a = zeros(N, 1);
m = ones(N, 3);

disp('Precomputations.');
tic;

% Each triangle is visited once.
for k = 1:N 
    nid = tri(k,:);
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    
    % Precompute phi's.
    lhs = [p(tri(k,:),:), ones(3,1)];
    phi(:,1,k) = lhs \ [1; 0; 0];
    phi(:,2,k) = lhs \ [0; 1; 0];
    phi(:,3,k) = lhs \ [0; 0; 1];
    
    % Precompute triangle area.
    J = [p1(1)-p3(1) p2(1)-p3(1);
         p1(2)-p3(2) p2(2)-p3(2)];
    a(k) = abs(det(J));
     
    % Precompute triangle midpoints.
    m(k,1:2) = 1/3 * (p1 + p2 + p3);
end
disp(['The loop took ' num2str(toc) 's.']);

% Initialize linear system.
M = sparse(n,n); % Mass matrix.
A_iso = sparse(n,n); % Isotropic stiffness matrix.
A_ani = sparse(n,n); % Anisotropic stiffness matrix.
b = zeros(n,1); % The right hand side.

disp('Started isotropic diffusion.');
tic;

% Each triangle is visited once.
for k = 1:N
    nid = tri(k,:);
    
    for i = 1:3
        for j = 1:3
            A_iso(nid(i),nid(j)) = A_iso(nid(i),nid(j)) + 0.5 * a(k) * (phi(1:2,i,k)' * phi(1:2,j,k));         
            M(nid(i),nid(j)) = M(nid(i),nid(j)) + 0.5 * a(k) * (phi(:,i,k)' * m(k,:)') * (phi(:,j,k)' * m(k,:)');
        end
    end
end
disp(['The loop took ' num2str(toc) 's.']);

u_iso = (M + sigma * A_iso) \ (M*u);

% Compute the diffusion coefficients.
d = zeros(N,2);

disp('Started computation of diffusion coefficients.');
tic;

% Each triangle is visited once.
for k = 1:N
    nid = tri(k,:);
    
    for i = 1:3
        d(k,:) = d(k,:) + u_iso(nid(i)) * phi(1:2,i,k)';
    end
end
disp(['The loop took ' num2str(toc) 's.']);

% Apply transfer function.
g_ani = g(sqrt(d(:,1).^2+d(:,2).^2));

disp('Started anisotropic diffusion.');
tic;

% Each triangle is visited once.
for k = 1:N
    nid = tri(k,:);
    for i = 1:3
        for j = 1:3
            A_ani(nid(i),nid(j)) = A_ani(nid(i),nid(j)) + 0.5 * a(k) * g_ani(k) * (phi(1:2,i)' * phi(1:2,j));
        end
    end
end
disp(['The loop took ' num2str(toc) 's.']);

u_ani = (M + tau * A_ani) \ (M*u);
