function [U, U_true] = loadCameraman(scale, sigma)

U = mat2gray(imread('cameraman.tif'));
U = imresize(U, scale);
U_true = U;

if sigma > 0
    U = U + sigma * randn(size(U));
    U = max(min(U,1),0);
end

disp(['Loaded camerman with ' num2str(size(U)) ' pixels.'])

end

