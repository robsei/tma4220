function imtriplot(tri, p, edge)

hold on;
triplot(tri, p(:,1), p(:,2));
plot(p(edge',1), p(edge',2), '-r', 'LineWidth', 2);

view(90,90)
axis off tight square
% title(['#nodes: ', num2str(size(p,1)), ' #elements: ', num2str(size(tri,1))]);

end

