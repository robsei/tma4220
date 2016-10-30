function iminterpsurf(tri, p, u, nx)

% Get the query points.
[X2q, X1q] = meshgrid(1:nx(1), 1:nx(2));

% Interpolate.
vq = griddata(p(:,1), p(:,2), u, X1q, X2q, 'natural');

% Plot.
imagesc(vq);
axis off tight square

end

