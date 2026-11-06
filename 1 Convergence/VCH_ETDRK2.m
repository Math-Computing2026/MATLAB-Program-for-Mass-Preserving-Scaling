












function [u] = VCH_ETDRK2(Nx,Ny,tau,T,kappa)


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
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + (abs(z)<=1e-14).*0.5;
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u0 = 0.5*sin(X).*sin(Y); 
u = u0; t = (0:Nt)'*tau;
for n = 1:Nt
    u_hat = fft2(u);
    
    N0 = f(u) - epsilon^2/nu*u; N0_hat = fft2(N0);
    Nu_hat = B_hat .* N0_hat + kappa * u_hat;
    Nu = real(ifft2(Nu_hat));

    w1_hat = phi0(z).*u_hat + tau*phi1(z).*Nu_hat;
    w1 = real(ifft2(w1_hat));
    N1 = f(w1) - epsilon^2/nu*w1; N1_hat = fft2(N1);
    Nw1_hat = B_hat .* N1_hat + kappa * w1_hat;
    Nw1 = real(ifft2(Nw1_hat));

    c1_1 = Nw1 - Nu; c1_1_hat = fft2(c1_1);
    u_new_hat = w1_hat + tau*phi2(z) .* c1_1_hat;
    u_new = real(ifft2(u_new_hat));
    u = u_new;
end
end

