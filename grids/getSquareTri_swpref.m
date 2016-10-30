function [p, tri, edge] = getSquareTri_swpref(n1, n2, u)

x1 = linspace(0,1,n1);
x2 = linspace(0,1,n2);

[X1_const, X2_const] = meshgrid(x1, x2);

p = [X1_const(:), X2_const(:)];
tri = [];

for i = 0:n1-2
    for j = 0:n2-2
        % Square corners:
        tl = (i)*n2 + (j) + 1;
        bl = (i+1)*n2 + (j) + 1;
        tr = (i)*n2 + (j+1) + 1;
        br = (i+1)*n2 + (j+1) + 1;
        
        if abs(u(tl) - u(br)) > abs(u(tr) - u(bl))
            tri(end+1,:) = [tl bl tr];
            tri(end+1,:) = [tr bl br];
        else
            tri(end+1,:) = [tl br tr];
            tri(end+1,:) = [tl bl br];
        end
    end
end

idx_diff = @(idx) [idx(1:end-1); idx(2:end)];

edge_n = idx_diff(1:n2)';
edge_s = idx_diff((n1-1)*n2+1:n1*n2)';
edge_w = idx_diff((0:n1-1)*n2+1)';
edge_e = idx_diff((1:n1)*n2)';

edge = [edge_n; edge_s; edge_w; edge_e];

end

