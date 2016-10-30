function I = quadrature1D(a,b,Nq,g)

% Function g should return a row vector.

if Nq == 1
    zq = 0;
    rhoq = 2;
elseif Nq == 2
    zq = [-sqrt(1/3) sqrt(1/3)];
    rhoq = [1 1];
elseif Nq == 3
    zq = [-sqrt(3/5) 0 sqrt(3/5)];
    rhoq = [5/9 8/9 5/9];
elseif Nq == 4
    zq = [-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];
    rhoq = [(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
else
    error('Number of quadrature points not supported');
end

% Change integration variables such that integral goes over [-1,1].
% x is the physical coordinates, z the reference coordinates.
x = @(z) (z+1)*(b-a)/2 + a;
h = @(z) g(x(z));

% Quadrature formula.
I = (b-a) / 2 * h(zq) * rhoq';

end

