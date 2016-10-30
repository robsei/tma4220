function plotdiffcoeff(nx1, nx2, g_ani)

odd = 1:2:length(g_ani);
even = odd + 1;

g_square = 0.5 * (g_ani(odd) + g_ani(even));

g_plotable = reshape(g_square, nx1-1, nx2-1);

imagesc(g_plotable);

end