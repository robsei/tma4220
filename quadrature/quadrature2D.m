function I = quadrature2D(p1,p2,p3,Nq,g)

% Vectors p1, p2, p3 have to be row vectors.
% Function g should accept one point per column and return a row vector.

if Nq == 1
    lbq = [1/3 1/3 1/3];
    rhoq = 1;
elseif Nq == 3
    lbq = [1/2 1/2 0;
          1/2 0 1/2;
          0 1/2 1/2];
    rhoq = [1/3 1/3 1/3];
elseif Nq == 4
    lbq = [1/3 1/3 1/3;
          3/5 1/5 1/5
          1/5 3/5 1/5
          1/5 1/5 3/5];
    rhoq = [-9/16 25/48 25/48 25/48];
else
    error('Number of quadrature points not supported');
end

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

end

