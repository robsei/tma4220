clear all; clc; close all;

% Load image.
U = mat2gray(imread('cameraman.tif'));
U = imresize(U, 0.3);

[n1, n2] = size(U);

[p tri edge] =  getSquareTri_swp(n1,n2,U(:));

edge_plotable = [edge ones(size(edge,1),1)];

subplot(1,2,1);
triplot(tri, p(:,1), p(:,2));
subplot(1,2,2);
triplot(edge_plotable, p(:,1), p(:,2));