function I = quadrature3D(p1,p2,p3,p4,Nq,g)

% Vectors p1, p2, p3, p4 have to be row vectors.
% Function g should accept one point per column and return a row vector.

if Nq == 1
    lbq = [1/4 1/4 1/4 1/4];
    rhoq = 1;
elseif Nq == 4
    lbq = [0.5854102 0.1381966 0.1381966 0.1381966;
          0.1381966 0.5854102 0.1381966 0.1381966;
          0.1381966 0.1381966 0.5854102 0.1381966;
          0.1381966 0.1381966 0.1381966 0.5854102];
    rhoq = [1/4 1/4 1/4 1/4];
elseif Nq == 5
    lbq = [1/4 1/4 1/4 1/4;
          1/2 1/6 1/6 1/6;
          1/6 1/2 1/6 1/6;
          1/6 1/6 1/2 1/6;
          1/6 1/6 1/6 1/2];
    rhoq = [-4/5 9/20 9/20 9/20 9/20];
else
    error('Number of quadrature points not supported');
end

% x is the physical coordinates written as a function of barycentric
% coordinates lb w.r.t to the (fixed) vertices.
% (Returns one 3d-point per row.)
x = @(lb) lb * [p1; p2; p3; p4];
% h is g in barycentric coordinates w.r.t. to the vertices.
h = @(lb) g(x(lb)');

% Jacobian of transformation function.
J = [p1(1)-p4(1) p2(1)-p4(1) p3(1)-p4(1);
    p1(2)-p4(2) p2(2)-p4(2) p3(2)-p4(2);
    p1(3)-p4(3) p2(3)-p4(3) p3(3)-p4(3)];

% Quadrature formula.
I = 1/6 * abs(det(J)) * h(lbq) * rhoq';

end

