















function [u] = VCH_ETD1(Nx,Ny,tau,T,kappa)


epsilon = 0.1; nu = 1; Nt = round(T/tau);
xa = 0; xb = 2*pi; hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; 
ya = 0; yb = 2*pi; hy = (yb-ya)/Ny; j = (0:Ny-1)'; y = ya + j*hy; 
[X,Y] = meshgrid(x,y);
f = @(u) u.^3 - u;

kx = [0:Nx/2-1,-Nx/2:-1]'; ky = [0:Ny/2-1,-Ny/2:-1]';
Lx = -4/hx^2 * sin(pi*kx/Nx).^2; Ly = -4/hy^2 * sin(pi*ky/Ny).^2;
[LX,LY] = meshgrid(Lx,Ly);
L_hat =  LX + LY;
Lk_hat = epsilon^2/nu * L_hat - kappa;

z = tau*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) ( exp(z)-1 ) ./ z;
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u0 = 0.5*sin(X).*sin(Y); 
u = u0; t = (0:Nt)'*tau;
for n = 1:Nt
    u_hat = fft2(u);
    
    N0 = f(u) - epsilon^2/nu*u;
    N0_hat = fft2(N0);
    Nu_hat = B_hat .* N0_hat + kappa * u_hat;
    
    u_new_hat = phi0(z) .* u_hat + tau * phi1(z) .* Nu_hat;
    u_new = real(ifft2(u_new_hat));
    u = u_new;
end

end

