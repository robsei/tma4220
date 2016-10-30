clear all; clc; close all;

% Load needed functions.
addpath('../quadrature');
addpath('../grids');

% Grid construction.
nx1 = 100;
nx2 = 100;
dx1 = 1/(nx1-1);
dx2 = 1/(nx2-1);
[p, tet, edge] = getUniform(nx1, nx2);

% Boundary conditions.
u = @(X) 1;

% Number of nodes.
n = nx1 * nx2;
% Number of elements.
N = size(tet, 1);

% Allocate memory.
A = zeros(n,n);
b = zeros(n,1);

% Each element is visited once.
for i=0:nx1-2
    for j=0:nx2-2
        A(i + j*nx1 + 1,i + j*nx1 + 1) = A(i + j*nx1 + 1,i + j*nx1 + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A(i + j*nx1 + 1,i + nx1*(j + 1) + 1) = A(i + j*nx1 + 1,i + nx1*(j + 1) + 1) + -dx1*dx2*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A(i + j*nx1 + 1,i + j*nx1 + 2) = A(i + j*nx1 + 1,i + j*nx1 + 2) + -dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A(i + j*nx1 + 1,i + nx1*(j + 1) + 2) = A(i + j*nx1 + 1,i + nx1*(j + 1) + 2) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + nx1*(j + 1) + 1,i + j*nx1 + 1) = A(i + nx1*(j + 1) + 1,i + j*nx1 + 1) + -dx1*dx2*(((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1) + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) = A(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 1) + dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)^2);
        A(i + nx1*(j + 1) + 1,i + j*nx1 + 2) = A(i + nx1*(j + 1) + 1,i + j*nx1 + 2) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) = A(i + nx1*(j + 1) + 1,i + nx1*(j + 1) + 2) + -dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A(i + j*nx1 + 2,i + j*nx1 + 1) = A(i + j*nx1 + 2,i + j*nx1 + 1) + -dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A(i + j*nx1 + 2,i + nx1*(j + 1) + 1) = A(i + j*nx1 + 2,i + nx1*(j + 1) + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + j*nx1 + 2,i + j*nx1 + 2) = A(i + j*nx1 + 2,i + j*nx1 + 2) + dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)^2);
        A(i + j*nx1 + 2,i + nx1*(j + 1) + 2) = A(i + j*nx1 + 2,i + nx1*(j + 1) + 2) + -dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + nx1*(j + 1) + 2,i + j*nx1 + 1) = A(i + nx1*(j + 1) + 2,i + j*nx1 + 1) + dx1*dx2*(((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2) + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) = A(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 1) + -dx1*dx2*((j/dx1 - (j + 1/2)/dx1)^2 + ((i + 1)/dx2 - (i + 1/2)/dx2)*(i/dx2 - (i + 1/2)/dx2));
        A(i + nx1*(j + 1) + 2,i + j*nx1 + 2) = A(i + nx1*(j + 1) + 2,i + j*nx1 + 2) + -dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + ((j + 1)/dx1 - (j + 1/2)/dx1)*(j/dx1 - (j + 1/2)/dx1));
        A(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) = A(i + nx1*(j + 1) + 2,i + nx1*(j + 1) + 2) + dx1*dx2*((i/dx2 - (i + 1/2)/dx2)^2 + (j/dx1 - (j + 1/2)/dx1)^2);
        b(i + j*nx1 + 1) = b(i + j*nx1 + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
        b(i + j*nx1 + 1) = b(i + j*nx1 + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + j*nx1 + 1) = b(i + j*nx1 + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        b(i + j*nx1 + 1) = b(i + j*nx1 + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 1) = b(i + nx1*(j + 1) + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 1) = b(i + nx1*(j + 1) + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))^2;
        b(i + nx1*(j + 1) + 1) = b(i + nx1*(j + 1) + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 1) = b(i + nx1*(j + 1) + 1) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + j*nx1 + 2) = b(i + j*nx1 + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        b(i + j*nx1 + 2) = b(i + j*nx1 + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + j*nx1 + 2) = b(i + j*nx1 + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
        b(i + j*nx1 + 2) = b(i + j*nx1 + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 2) = b(i + nx1*(j + 1) + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 2) = b(i + nx1*(j + 1) + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 2) = b(i + nx1*(j + 1) + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * -dx1*dx2*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        b(i + nx1*(j + 1) + 2) = b(i + nx1*(j + 1) + 2) + u([(i+1/2)*dx1, (j+1/2)*dx2]) * dx1*dx2*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2))^2;
    end
end

% Treat boundary conditions.
A(edge(:,1),:) = [];
A(:,edge(:,1)) = [];
b(edge(:,1)) = [];

% Solve the system.
u = A \ b;

% Fill in boundary condition.
U = zeros(n,1);
idx = 1:n;
idx(edge(:,1)) = [];
U(idx) = u;

% Plot the numerical result.
imagesc(reshape(U,nx1,nx2));