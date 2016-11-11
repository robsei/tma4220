clear all; clc; close all;

p1 = [0 0];
p2 = [1 3];
Nqs = 1:4;
g = @(x) x(1,:).^2 + x(2,:).^4;

I_true = sqrt(10) * (1/3 + (3^4)/5);

for Nq = Nqs
    I = linequadrature1D(p1,p2,Nq,g);
    disp(['For Nq=' num2str(Nq) ': I = ' num2str(I) ', |I - I_true| = ' num2str(abs(I - I_true))]);
end