

function [Mass_Err] = Mass_ETD2(N,xa,xb,tau,T,epsilon,u0)


nu = 1; Nt = round(T/tau); Nx = N; Ny = N;
hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; hy = hx; y = x;
[X,Y] = meshgrid(x,y);

theta = 0.8; theta_c = 1.6; u_safe = @(u) max(min(u, 1-1e-6), -1+1e-6);
f = @(u) theta/2*log((1+u_safe(u))./(1-u_safe(u))) - theta_c*u; 
beta = 0.98678947;  kappa = 28.88; 

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

u = u0; t = (0:Nt)'*tau;
Mass = zeros(Nt+1,1); Mass_Err = zeros(Nt,1);  
Mass(1) = hx*hy*sum(u0(:));
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
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
    Mass(n+1,1) = hx*hy*sum(u(:)); Mass_Err(n,1) = abs(Mass(n+1,1)-Mass(1));
end
end

