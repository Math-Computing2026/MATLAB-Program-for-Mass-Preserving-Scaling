
function [u, u_hist, t] = MBP_ETDRK2_3D(N, xa, xb, tau, T, epsilon, u0)

nu = 1; Nt = round(T/tau); Nx = N; Ny = N; Nz = N;
hx = (xb-xa)/Nx; hy = hx; hz = hx; 
f = @(u) u.^3 - u; beta = 1; kappa = 3;




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


u = u0; t = (0:Nt)'*tau;
u_hist = zeros(Nt+1, 1); u_hist(1) = max(abs(u0(:)));
for n = 1:Nt
    u_hat = fftn(u);
    
    F0 = f(u) - epsilon^2/nu*u; 
    F0_hat = fftn(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifftn(Nu_hat));

    Nu = max(min(Nu, kappa*beta), -kappa*beta); 
    Nu_hat = fftn(Nu);
    

    w1_hat = phi0(z).*u_hat + tau*phi1(z).*Nu_hat;
    w1 = real(ifftn(w1_hat));
    
    N1 = f(w1) - epsilon^2/nu*w1; 
    N1_hat = fftn(N1);
    Nw1_hat = B_hat .* N1_hat + kappa * w1_hat;
    

    c1_1_hat = Nw1_hat - Nu_hat; 
    u_new_hat = w1_hat + tau*phi2(z) .* c1_1_hat;
    u_new = real(ifftn(u_new_hat));
    u = u_new;   u_hist(n+1, 1) = max(abs(u(:)));
  
end
end