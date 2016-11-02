function [p, tri, edge] = getSquareTri_swp(U)

u = U(:);
nx = size(U);

x1 = 1:nx(1);
x2 = 1:nx(2);

[X1_const, X2_const] = meshgrid(x1, x2);

p = [X1_const(:), X2_const(:)];
tri = [];

for i = 0:nx(1)-2
    for j = 0:nx(2)-2
        % Square's corners:
        tl = (i)*nx(2) + (j) + 1;
        bl = (i+1)*nx(2) + (j) + 1;
        tr = (i)*nx(2) + (j+1) + 1;
        br = (i+1)*nx(2) + (j+1) + 1;
        
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

edge_n = idx_diff(1:nx(2))';
edge_s = idx_diff((nx(1)-1)*nx(2)+1:nx(1)*nx(2))';
edge_w = idx_diff((nx(1)-1:-1:0)*nx(2)+1)';
edge_e = idx_diff((1:nx(1))*nx(2))';

edge = [edge_n; edge_e; edge_s; edge_w];

end

