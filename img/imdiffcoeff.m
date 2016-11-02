function imdiffcoeff(g_ani, nx)
% This function visualizes the g_ani on a *structured grid only*.

% The odd index gets all triangles with top left corner as member.
odd = 1:2:length(g_ani);
% The even index gets all triangles with bottom right corner as member.
even = odd + 1;

% Take the average of the two triangles.
g_square = 0.5 * (g_ani(odd) + g_ani(even));

% Go back to rectangular image format.
G = reshape(g_square, nx(1)-1, nx(2)-1);

% Visualize.
imagesc(G);
view(-90,-90)
axis off tight square

end