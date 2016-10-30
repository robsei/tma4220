clear all; clc; close all;

addpath('../img');

% Load image.
U = loadCameraman(1, 0.01);

[p tri edge p_idx] =  getVarigrid(U, 0.001, 50);

u_plottable = U(sub2ind(size(U), p_idx(:,1), p_idx(:,2)));

edge_plotable = [edge ones(size(edge,1),1)];

subplot(2,2,1);
imshow(U);
subplot(2,2,2);
triplot(tri, p(:,1), p(:,2));
subplot(2,2,3);
imtrisurf(tri, p, u_plottable);
subplot(2,2,4);
triplot(edge_plotable, p(:,1), p(:,2));