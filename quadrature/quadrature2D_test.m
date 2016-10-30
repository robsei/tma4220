% clear all; clc; close all;

p1 = [1 0];
p2 = [3 1];
p3 = [3 2];
Nqs = [1 3 4];
g = @(x) log(x(1,:) + x(2,:));
% g = @(x) x(1,:) + x(2,:);
I_true = 1.16542;
% I_true = 10/3;

for Nq = Nqs
    I = quadrature2D(p1,p2,p3,Nq,g);
    disp(['For Nq=' num2str(Nq) ': I = ' num2str(I) ', |I - I_true| = ' num2str(abs(I - I_true))]);
end

% Command for www.wolframalpha.com:
% integrate log(x+y) dy dx x=1..3 y=1/2(x-1)..x-1
