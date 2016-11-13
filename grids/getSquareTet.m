function [p tet edge] = getSquareTet(nx1, nx2)

disp('Started grid generation.');
tic;

x1 = linspace(0,1,nx1);
x2 = linspace(0,1,nx2);

sub2ind = @(i,j) i + nx1 * j + 1;

[X1_const, X2_const] = meshgrid(x1, x2);

p = [X1_const(:), X2_const(:)];

tet = [];
for i = 0:nx1-2
    for j = 0:nx2-2
        tet(end+1,:) = [sub2ind(i,j) sub2ind(i+1,j) sub2ind(i,j+1) sub2ind(i+1,j+1)];
    end
end

edge = [];
for i = 0:nx1-1
    j = 0;
    edge(end+1,:) = [sub2ind(i,j) sub2ind(i+1,j)];
    j = nx2-1;
    edge(end+1,:) = [sub2ind(i,j) sub2ind(i+1,j)];
end
for j = 0:nx2-1
    i = 0;
    edge(end+1,:) = [sub2ind(i,j) sub2ind(i,j+1)];
    i = nx1-1;
    edge(end+1,:) = [sub2ind(i,j) sub2ind(i,j+1)];
end

disp(['Grid generation took ' num2str(toc) 's.']);
disp(['Number of nodes: ' num2str(size(p,1))])
disp(['Number of elements: ' num2str(size(tet,1))])

end

