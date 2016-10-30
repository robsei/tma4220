function U = loadCameraman(scale, sigma)

U = mat2gray(imread('cameraman.tif'));
U = imresize(U, scale);

if sigma > 0
    U = U + sigma * randn(size(U));
    U = max(min(U,1),0);
end

end

