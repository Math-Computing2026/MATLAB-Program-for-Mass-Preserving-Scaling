
function [u, E, t] = Energy_ETDRK2_3D(N, xa, xb, tau, T, epsilon, u0)

Nx = N; Ny = N; Nz = N; nu = 1; Nt = round(T/tau); 
hx = (xb-xa)/Nx;  hy = hx; hz = hx; 

F = @(u) 1/4*(1-u.^2).^2; f = @(u) u.^3 - u; 

theta = 0.8; theta_c = 1.6; u_safe = @(u) max(min(u, 1-1e-6), -1+1e-6);
F = @(u) theta/2*( (1+u_safe(u)).*log(1+u_safe(u)) + (1-u_safe(u)).*log(1-u_safe(u)) ) - theta_c/2*u.^2;
f = @(u) theta/2*log((1+u_safe(u))./(1-u_safe(u))) - theta_c*u; 
beta = 0.986807;  kappa = 10;


kx = [0:Nx/2-1, -Nx/2:-1]'; ky = [0:Ny/2-1, -Ny/2:-1]';kz = [0:Nz/2-1, -Nz/2:-1]';
Lx = -4/hx^2 * sin(pi*kx/Nx).^2; Ly = -4/hy^2 * sin(pi*ky/Ny).^2; Lz = -4/hz^2 * sin(pi*kz/Nz).^2;
[LX, LY, LZ] = ndgrid(Lx, Ly, Lz);
L_hat = LX + LY + LZ;  
Lk_hat = epsilon^2/nu * L_hat - kappa;


z = tau*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + (abs(z)<=1e-14).*0.5;
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);


u = u0; t = (0:Nt)'*tau; E = zeros(Nt+1, 1); 
E(1) = Energy_FD_3D(u0, epsilon, F, hx, hy, hz);
for n = 1:Nt
    u_hat = fftn(u);
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fftn(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifftn(Nu_hat));

    Nu = max(min(Nu, kappa*beta), -kappa*beta); Nu_hat = fftn(Nu);
    

    w1_hat = phi0(z).*u_hat + tau*phi1(z).*Nu_hat;
    w1 = real(ifftn(w1_hat));
    N1 = f(w1) - epsilon^2/nu*w1; 
    N1_hat = fftn(N1);
    Nw1_hat = B_hat .* N1_hat + kappa * w1_hat;
    


    c1_1_hat = Nw1_hat - Nu_hat;   
    u_new_hat = w1_hat + tau*phi2(z) .* c1_1_hat;
    u_new = real(ifftn(u_new_hat));
    u = u_new;
  

    E(n+1, 1) = Energy_FD_3D(u, epsilon, F, hx, hy, hz);
end
end

function E_FD = Energy_FD_3D(u, epsilon, F, hx, hy, hz)

ux = (circshift(u, [-1, 0, 0]) - circshift(u, [1, 0, 0])) / (2*hx);
uy = (circshift(u, [0, -1, 0]) - circshift(u, [0, 1, 0])) / (2*hy);
uz = (circshift(u, [0, 0, -1]) - circshift(u, [0, 0, 1])) / (2*hz);

energy_density = 1/2 * epsilon^2 * (ux.^2 + uy.^2 + uz.^2) + F(u);
E_FD = hx * hy * hz * sum(energy_density(:));
end