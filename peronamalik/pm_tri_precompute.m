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

lbq = [1/3 1/3 1/3];
rhoq = 1;

% x is the physical coordinates written as a function of barycentric
% coordinates lb w.r.t to the (fixed) vertices.
% (Returns one 2d-point per row.)
x = @(lb) lb * [p1; p2; p3];
% h is g in barycentric coordinates w.r.t. to the vertices.
h = @(lb) g(x(lb)');

% Jacobian of transformation function.
J = [p1(1)-p3(1) p2(1)-p3(1);
     p1(2)-p3(2) p2(2)-p3(2)];

% Quadrature formula.
I = 0.5 * abs(det(J)) * h(lbq) * rhoq';
