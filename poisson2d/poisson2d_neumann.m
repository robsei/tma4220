function u_bd = poisson2d_neumann(p, tri, edge, f, g, eps)

% Load needed functions.
addpath('../quadrature');

% Number of points.
n = size(p, 1);

% Number of elements.
N = size(tri, 1);

% Prepare linear system.
A = sparse(n,n);
b = zeros(n,1);

% Order of integration.
Nq = 1;

% Get (Dirichlet) boundary nodes.
e = unique(edge(:));

% Visit each triangle once.
for k = 1:N
    % Get node IDs to identify points and basis functions.
    nid = tri(k,:);
    p1 = p(nid(1),:);
    p2 = p(nid(2),:);
    p3 = p(nid(3),:);
    
    % Compute PHIs that have support in this triangle.
    % Create LHS of the linear system.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    % Integrate and update A and b.
    for i = 1:3
        for j = 1:3
            g1 = @(x) dot(phi(1:2,i), phi(1:2,j)) * ones(1, size(x,2));
            A(nid(i),nid(j)) = A(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,g1);
        end
        
        g2 = @(x) f(x) * phi(:,i)' * [x; ones(1, size(x,2))];
        b(nid(i)) = b(nid(i)) + quadrature2D(p1,p2,p3,Nq,g2);
    end
    
    % Treatment of Neumann boundary conditions.
    % Which points are on Omega_N?
    is_bd_point = ...
        (abs(p(nid,1)) < eps & p(nid,2) < eps & p(nid,2) > -1) ...
             | ...
        (p(nid,1) > eps & abs(p(nid,2)) < eps & p(nid,1) < 1);
    
    % If two points are on Omage_N, triangle has an edge with Neumann boundary.
    has_nm_edge = sum(is_bd_point) == 2;
    
    % Update e and b, if necessary.
    if has_nm_edge
        nm_edge_nid = nid(is_bd_point);
        nm_edge_pt = p(nm_edge_nid,:);
    
        for i = find(is_bd_point')
            e(e == nid(i)) = [];
            g3 = @(X) g(X) * phi(:,i)' * [X; 1];
            b(nid(i)) = b(nid(i)) + linequadrature1D(nm_edge_pt(1,:),nm_edge_pt(2,:), Nq, g3);
        end
    end
end

% Get indices of inner nodes and Neumann nodes.
idx = 1:n;
idx(e) = [];

% Solve linear system on inner nodes and Neumann nodes.
u = A(idx,idx) \ b(idx);

% Fill in Dirichlet boundary condition.
u_bd = zeros(n,1);
u_bd(idx) = u;