function [U, U_true]  = loadRhombus(scale, sigma)

nx = floor(scale/4)*2*2;

U = zeros(nx);

pos_rhombus = nx * [1/4 2/4 2/4 1/4 3/4 2/4 2/4 3/4];
U = insertShape(U,'FilledPolygon', pos_rhombus, 'Color', 'white', 'Opacity', 1, 'SmoothEdges', false);
U = U(:,:,1);

U_true = U;

if sigma > 0
    U = U + sigma * randn(size(U));
    U = max(min(U,1),0);
end

disp(['Loaded rhombus with ' num2str(size(U)) ' pixels.'])

end

