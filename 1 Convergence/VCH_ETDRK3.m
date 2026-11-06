


clear all; clc;
Nx = 128; Ny = Nx; T = 1; kappa = 3;
for i = 1:7
   tau = 2^(-i);
   [u] = VCH_ETDRK3_2D(Nx,Ny,tau,T,kappa);
   [U] = VCH_ETDRK3_2D(Nx,Ny,tau/2,T,kappa);
   MaxErr(i) = max(max(abs(U-u)));
end
Ord_time = log2(MaxErr(1:end-1)./MaxErr(2:end))

function [u] = VCH_ETDRK3_2D(Nx,Ny,tau,T,kappa)


epsilon = 0.01; nu = 1; Nt = round(T/tau);
xa = 0; xb = 2*pi; hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; 
ya = 0; yb = 2*pi; hy = (yb-ya)/Ny; j = (0:Ny-1)'; y = ya + j*hy; 
[X,Y] = meshgrid(x,y);
f = @(u) u.^3 - u;

kx = [0:Nx/2-1,-Nx/2:-1]'; ky = [0:Ny/2-1,-Ny/2:-1]';
Lx = -4/hx^2 * sin(pi*kx/Nx).^2; Ly = -4/hy^2 * sin(pi*ky/Ny).^2; 
[LX,LY] = meshgrid(Lx,Ly);
L_hat =  LX + LY;
Lk_hat = epsilon^2/nu * L_hat - kappa;

z = tau*Lk_hat; z_half = tau/2*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + 1/2.*(abs(z)<=1e-14);
phi3 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2)./(z.^3)) + 1/6.*(abs(z)<=1e-14);
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u0 = 0.5*sin(X).*sin(Y);
u = u0; t = (0:Nt)'*tau;
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifft2(Nu_hat));

    w1_hat = phi0(z).*u_hat + tau*phi1(z).*Nu_hat;
    w1_half_hat = phi0(z_half).*u_hat + tau/2*phi1(z_half).*Nu_hat;
    w1 = real(ifft2(w1_hat));
    F1 = f(w1) - epsilon^2/nu*w1; F1_hat = fft2(F1);
    Nw1_hat = B_hat .* F1_hat + kappa * w1_hat;
    Nw1 = real(ifft2(Nw1_hat));

    c1_1 = Nw1 - Nu; c1_1_hat = fft2(c1_1);
    w2_hat = w1_hat + tau*phi2(z) .* c1_1_hat; 
    w2 = real(ifft2( w2_hat));
    w2_half_hat = w1_half_hat + tau*(1/2)^2*phi2(z_half).*c1_1_hat;
    w2_half = real(ifft2( w2_half_hat));
    
    F2_half = f(w2_half) - epsilon^2/nu*w2_half; F2_half_hat = fft2(F2_half);
    F2 = f(w2) - epsilon^2/nu*w2; F2_hat = fft2(F2);
    
    Nw2_s1_hat = B_hat .* F2_half_hat + kappa * w2_half_hat;
    Nw2_s2_hat = B_hat .* F2_hat + kappa * w2_hat;
    Nw2_s1 = real(ifft2(Nw2_s1_hat));
    Nw2_s2 = real(ifft2(Nw2_s2_hat));
    

    b1 = Nw2_s1 - Nu; b2 = Nw2_s2 - Nu;
    c2_1 = 4*b1 - b2; c2_2 = -4*b1 + 2*b2;
    c2_1_hat = fft2(c2_1); c2_2_hat = fft2(c2_2);
    
    u_new_hat = w1_hat + tau*phi2(z).*c2_1_hat + 2*tau*phi3(z).*c2_2_hat;
    u_new = real(ifft2(u_new_hat));
    u = u_new;
end
end

