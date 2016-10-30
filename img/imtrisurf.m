function imtrisurf(tri, p, u)

trisurf(tri, p(:,1), p(:,2), u, 'edgecolor','none');
colormap gray
view(90,90)
shading interp
axis square

end

