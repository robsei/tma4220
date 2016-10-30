clear all; clc; close all;

a = 1;
b = 2;
Nqs = 1:4;
g = @(x) exp(x);
I_true = exp(2) - exp(1);

for Nq = Nqs
    I = quadrature1D(a,b,Nq,g);
    disp(['For Nq=' num2str(Nq) ': I = ' num2str(I) ', |I - I_true| = ' num2str(abs(I - I_true))]);
end