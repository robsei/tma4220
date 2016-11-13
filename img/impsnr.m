function z = impsnr(p, u_ani, U_true)

nx = size(U_true);

% Get the query points.
[X2q, X1q] = meshgrid(1:nx(1), 1:nx(2));

% Interpolate.
vq = griddata(p(:,1), p(:,2), u_ani, X1q, X2q, 'natural');

z = psnr(vq, U_true);

end

