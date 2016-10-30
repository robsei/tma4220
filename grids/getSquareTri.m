function [p, tri, edge] = getSquareTri(n1, n2)

x1 = linspace(0,1,n1);
x2 = linspace(0,1,n2);

[X1_const, X2_const] = meshgrid(x1, x2);

p = [X1_const(:), X2_const(:)];
tri = [];

for i = 0:n1-2
    for j = 0:n2-2
        tri(end+1,:) = [i*n2 (i+1)*n2 i*n2+1] + j + 1;
        tri(end+1,:) = [(i+1)*n2 i*n2+1 (i+1)*n2+1] + j + 1;
    end
end

idx_diff = @(idx) [idx(1:end-1); idx(2:end)];

edge_n = idx_diff(1:n2)';
edge_s = idx_diff((n1-1)*n2+1:n1*n2)';
edge_w = idx_diff((0:n1-1)*n2+1)';
edge_e = idx_diff((1:n1)*n2)';

edge = [edge_n; edge_s; edge_w; edge_e];

end

