function [p, tri, edge] = getSquareTri(nx)

x1 = 1:nx(1);
x2 = 1:nx(2);

[X1_const, X2_const] = meshgrid(x1, x2);

p = [X1_const(:), X2_const(:)];
tri = [];

for i = 0:nx(1)-2
    for j = 0:nx(2)-2
        tri(end+1,:) = [i*nx(2) (i+1)*nx(2) i*nx(2)+1] + j + 1;
        tri(end+1,:) = [(i+1)*nx(2) i*nx(2)+1 (i+1)*nx(2)+1] + j + 1;
    end
end

idx_diff = @(idx) [idx(1:end-1); idx(2:end)];

edge_n = idx_diff(1:nx(2))';
edge_s = idx_diff((nx(1)-1)*nx(2)+1:nx(1)*nx(2))';
edge_w = idx_diff((0:nx(1)-1)*nx(2)+1)';
edge_e = idx_diff((1:nx(1))*nx(2))';

edge = [edge_n; edge_s; edge_w; edge_e];

end
