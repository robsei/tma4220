function [U_iso_bd, G_ani, U_ani_bd] = pm_tet_pre(U, sigma, tau, g)

% Grid generation.
[nx1, nx2] = size(U);
[~, ~, edge] = getSquareTet(nx1, nx2);

% Number of nodes.
n = nx1 * nx2;

% Allocate memory.
M = sparse(n,n); % Mass matrix.
A_iso = sparse(n,n); % Isotropic stiffness matrix.
A_ani = sparse(n,n); % Anisotropic stiffness matrix.
b = zeros(n,1); % The right hand side.
d = zeros(n, 2); % The gradients.

disp('Started isotropic diffusion.');
tic;

% Loop for isotropic diffusion.
% Each square is visited once.
for i = 1:nx1-1
    for j = 1:nx2-1
        A_iso(i + nx1*(j - 1),i + nx1*(j - 1)) = A_iso(i + nx1*(j - 1),i + nx1*(j - 1)) + 1/2;
        A_iso(i + nx1*(j - 1),i + j*nx1) = A_iso(i + nx1*(j - 1),i + j*nx1) + 0;
        A_iso(i + nx1*(j - 1),i + nx1*(j - 1) + 1) = A_iso(i + nx1*(j - 1),i + nx1*(j - 1) + 1) + 0;
        A_iso(i + nx1*(j - 1),i + j*nx1 + 1) = A_iso(i + nx1*(j - 1),i + j*nx1 + 1) + -1/2;
        A_iso(i + j*nx1,i + nx1*(j - 1)) = A_iso(i + j*nx1,i + nx1*(j - 1)) + 0;
        A_iso(i + j*nx1,i + j*nx1) = A_iso(i + j*nx1,i + j*nx1) + 1/2;
        A_iso(i + j*nx1,i + nx1*(j - 1) + 1) = A_iso(i + j*nx1,i + nx1*(j - 1) + 1) + -1/2;
        A_iso(i + j*nx1,i + j*nx1 + 1) = A_iso(i + j*nx1,i + j*nx1 + 1) + 0;
        A_iso(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) = A_iso(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) + 0;
        A_iso(i + nx1*(j - 1) + 1,i + j*nx1) = A_iso(i + nx1*(j - 1) + 1,i + j*nx1) + -1/2;
        A_iso(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) = A_iso(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) + 1/2;
        A_iso(i + nx1*(j - 1) + 1,i + j*nx1 + 1) = A_iso(i + nx1*(j - 1) + 1,i + j*nx1 + 1) + 0;
        A_iso(i + j*nx1 + 1,i + nx1*(j - 1)) = A_iso(i + j*nx1 + 1,i + nx1*(j - 1)) + -1/2;
        A_iso(i + j*nx1 + 1,i + j*nx1) = A_iso(i + j*nx1 + 1,i + j*nx1) + 0;
        A_iso(i + j*nx1 + 1,i + nx1*(j - 1) + 1) = A_iso(i + j*nx1 + 1,i + nx1*(j - 1) + 1) + 0;
        A_iso(i + j*nx1 + 1,i + j*nx1 + 1) = A_iso(i + j*nx1 + 1,i + j*nx1 + 1) + 1/2;
        M(i + nx1*(j - 1),i + nx1*(j - 1)) = M(i + nx1*(j - 1),i + nx1*(j - 1)) + ((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
        M(i + nx1*(j - 1),i + j*nx1) = M(i + nx1*(j - 1),i + j*nx1) + -((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j - 1),i + nx1*(j - 1) + 1) = M(i + nx1*(j - 1),i + nx1*(j - 1) + 1) + -((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j - 1),i + j*nx1 + 1) = M(i + nx1*(j - 1),i + j*nx1 + 1) + ((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1,i + nx1*(j - 1)) = M(i + j*nx1,i + nx1*(j - 1)) + -((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1,i + j*nx1) = M(i + j*nx1,i + j*nx1) + (j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))^2;
        M(i + j*nx1,i + nx1*(j - 1) + 1) = M(i + j*nx1,i + nx1*(j - 1) + 1) + (i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1,i + j*nx1 + 1) = M(i + j*nx1,i + j*nx1 + 1) + -(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) = M(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) + -((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j - 1) + 1,i + j*nx1) = M(i + nx1*(j - 1) + 1,i + j*nx1) + (i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) = M(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) + (i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))^2;
        M(i + nx1*(j - 1) + 1,i + j*nx1 + 1) = M(i + nx1*(j - 1) + 1,i + j*nx1 + 1) + -(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + nx1*(j - 1)) = M(i + j*nx1 + 1,i + nx1*(j - 1)) + ((i + 1)*(j + 1) - (i + 1)*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + j*nx1) = M(i + j*nx1 + 1,i + j*nx1) + -(j*(i + 1) - j*(i + 1/2) - (i + 1)*(j + 1/2) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + nx1*(j - 1) + 1) = M(i + j*nx1 + 1,i + nx1*(j - 1) + 1) + -(i*(j + 1) - i*(j + 1/2) - (i + 1/2)*(j + 1) + (i + 1/2)*(j + 1/2))*(i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2));
        M(i + j*nx1 + 1,i + j*nx1 + 1) = M(i + j*nx1 + 1,i + j*nx1 + 1) + (i*j - i*(j + 1/2) - j*(i + 1/2) + (i + 1/2)*(j + 1/2))^2;
    end
end
disp(['The loop took ' num2str(toc) 's.']);

% Compute RHS only once.
b = M*U(:);

% Solve linear system.
u_iso = (M + sigma * A_iso) \ b;
u_iso_bd = u_iso;

disp('Started computation of diffusion coefficients.');
tic;

for i = 1:nx1-1
    for j = 1:nx2-1
        d(i + nx1*(j - 1),1) = d(i + nx1*(j - 1),1) + -u_iso_bd(i + nx1*(j - 1))/2;
        d(i + nx1*(j - 1),2) = d(i + nx1*(j - 1),2) + -u_iso_bd(i + nx1*(j - 1))/2;
        d(i + nx1*(j - 1),1) = d(i + nx1*(j - 1),1) + -u_iso_bd(i + j*nx1)/2;
        d(i + nx1*(j - 1),2) = d(i + nx1*(j - 1),2) + u_iso_bd(i + j*nx1)/2;
        d(i + nx1*(j - 1),1) = d(i + nx1*(j - 1),1) + u_iso_bd(i + nx1*(j - 1) + 1)/2;
        d(i + nx1*(j - 1),2) = d(i + nx1*(j - 1),2) + -u_iso_bd(i + nx1*(j - 1) + 1)/2;
        d(i + nx1*(j - 1),1) = d(i + nx1*(j - 1),1) + u_iso_bd(i + j*nx1 + 1)/2;
        d(i + nx1*(j - 1),2) = d(i + nx1*(j - 1),2) + u_iso_bd(i + j*nx1 + 1)/2;
    end
end
disp(['The loop took ' num2str(toc) 's.']);

g_ani = g(sqrt(d(:,1).^2+d(:,2).^2));

disp('Started anisotropic diffusion.');
tic;

for i = 1:nx1-1
    for j = 1:nx2-1
        A_ani(i + nx1*(j - 1),i + nx1*(j - 1)) = A_ani(i + nx1*(j - 1),i + nx1*(j - 1)) + g_ani(i + nx1*(j - 1))/2;
        A_ani(i + nx1*(j - 1),i + j*nx1) = A_ani(i + nx1*(j - 1),i + j*nx1) + 0;
        A_ani(i + nx1*(j - 1),i + nx1*(j - 1) + 1) = A_ani(i + nx1*(j - 1),i + nx1*(j - 1) + 1) + 0;
        A_ani(i + nx1*(j - 1),i + j*nx1 + 1) = A_ani(i + nx1*(j - 1),i + j*nx1 + 1) + -g_ani(i + nx1*(j - 1))/2;
        A_ani(i + j*nx1,i + nx1*(j - 1)) = A_ani(i + j*nx1,i + nx1*(j - 1)) + 0;
        A_ani(i + j*nx1,i + j*nx1) = A_ani(i + j*nx1,i + j*nx1) + g_ani(i + nx1*(j - 1))/2;
        A_ani(i + j*nx1,i + nx1*(j - 1) + 1) = A_ani(i + j*nx1,i + nx1*(j - 1) + 1) + -g_ani(i + nx1*(j - 1))/2;
        A_ani(i + j*nx1,i + j*nx1 + 1) = A_ani(i + j*nx1,i + j*nx1 + 1) + 0;
        A_ani(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) = A_ani(i + nx1*(j - 1) + 1,i + nx1*(j - 1)) + 0;
        A_ani(i + nx1*(j - 1) + 1,i + j*nx1) = A_ani(i + nx1*(j - 1) + 1,i + j*nx1) + -g_ani(i + nx1*(j - 1))/2;
        A_ani(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) = A_ani(i + nx1*(j - 1) + 1,i + nx1*(j - 1) + 1) + g_ani(i + nx1*(j - 1))/2;
        A_ani(i + nx1*(j - 1) + 1,i + j*nx1 + 1) = A_ani(i + nx1*(j - 1) + 1,i + j*nx1 + 1) + 0;
        A_ani(i + j*nx1 + 1,i + nx1*(j - 1)) = A_ani(i + j*nx1 + 1,i + nx1*(j - 1)) + -g_ani(i + nx1*(j - 1))/2;
        A_ani(i + j*nx1 + 1,i + j*nx1) = A_ani(i + j*nx1 + 1,i + j*nx1) + 0;
        A_ani(i + j*nx1 + 1,i + nx1*(j - 1) + 1) = A_ani(i + j*nx1 + 1,i + nx1*(j - 1) + 1) + 0;
        A_ani(i + j*nx1 + 1,i + j*nx1 + 1) = A_ani(i + j*nx1 + 1,i + j*nx1 + 1) + g_ani(i + nx1*(j - 1))/2;
    end
end
disp(['The loop took ' num2str(toc) 's.']);

% Solve linear system.
u_ani = (M + tau * A_ani) \ b;
u_ani_bd = u_ani;

% Reshape all return values.
U_iso_bd = reshape(u_iso_bd, nx1, nx2);
G_ani = reshape(g_ani, nx1, nx2);
U_ani_bd = reshape(u_ani_bd, nx1, nx2);