function u_bd = poisson2d_dirichlet(p, tri, edge, f, Nq)

% Load needed functions.
addpath('../quadrature');

% Number of points.
n = size(p, 1);

% Number of elements.
N = size(tri, 1);

% Prepare linear system.
A = sparse(n,n);
b = zeros(n,1);

% Visit each triangle once.
for k = 1:N
    % Get node IDs to identify points and basis functions.
    nid = tri(k,:);
    
    % Compute PHIs that have support in this triangle.
    % Create LHS of the linear system.
    lhs = [p(nid,:), ones(3,1)];
    phi = zeros(3,3);
    phi(:,1) = lhs \ [1; 0; 0];
    phi(:,2) = lhs \ [0; 1; 0];
    phi(:,3) = lhs \ [0; 0; 1];
    
    % Integrate and update A and b.
    for i = 1:3
        p1 = p(nid(1),:);
        p2 = p(nid(2),:);
        p3 = p(nid(3),:);
        for j = 1:3
            g1 = @(x) dot(phi(1:2,i), phi(1:2,j)) * ones(1, size(x,2));
            A(nid(i),nid(j)) = A(nid(i),nid(j)) + quadrature2D(p1,p2,p3,Nq,g1);
        end
        g2 = @(x) f(x) * phi(:,i)' * [x; ones(1, size(x,2))];
        b(nid(i)) = b(nid(i)) +  quadrature2D(p1,p2,p3,Nq,g2);
    end
end

% Get boundary nodes.
e = unique(edge(:));

% Get indices of inner points.
idx = 1:n;
idx(e) = [];

% Solve linear system on inner points.
u = A(idx,idx) \ b(idx);

% Fill in boundary condition.
u_bd = zeros(n,1);
u_bd(idx) = u;

end