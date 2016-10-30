clear all; clc; close all;

p1 = [0 0 0];
p2 = [0 2 0];
p3 = [0 0 2];
p4 = [2 0 0];
Nqs = [1 4 5];
g = @(x) exp(x(1,:));
I_true = 2 + 4 + 4 * (exp(2) - 1);
I_true = 2 * exp(2) - 10;
I_true = 1.79179;

% Other example:
% http://mathinsight.org/triple_integral_examples

p1 = [0 0 0];
p2 = [2 0 0];
p3 = [0 3 0];
p4 = [0 0 1];
Nqs = [1 4 5];
g = @(x) ones(1,size(x,2));
I_true = 1;
g = @(x) exp(x(1,:));
I_true = 1.79179;

for Nq = Nqs
    I = quadrature3D(p1,p2,p3,p4,Nq,g);
    disp(['For Nq=' num2str(Nq) ': I = ' num2str(I) ', |I - I_true| = ' num2str(abs(I - I_true))]);
end