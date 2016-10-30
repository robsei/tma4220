clear all; clc; close all;

addpath('../img');

% Load image.
U = loadCameraman(1, 0.05);

[p, tri, edge] =  getGradgrid(U, 2, @(x) power(x,1));
u = maskApply(U,p);

subplot(2,2,1);
imshow(U);
subplot(2,2,2);
imtriplot(tri, p, edge);
subplot(2,2,3);
imtrisurf(tri, p, u);
subplot(2,2,4);
iminterpsurf(tri, p, u, size(U));

print('../out/test_getGradgrid','-dpng');

