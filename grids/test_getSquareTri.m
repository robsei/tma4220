clear all; clc; close all;

n1 = 6;
n2 = 8;

n1 = 77;
n2 = 77;

[p tri edge] =  getSquareTri(n1,n2);

edge_plotable = [edge ones(size(edge,1),1)];

subplot(1,2,1);
triplot(tri, p(:,1), p(:,2));
subplot(1,2,2);
triplot(edge_plotable, p(:,1), p(:,2));