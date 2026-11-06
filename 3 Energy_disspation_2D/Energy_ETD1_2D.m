


function [u,E,t] = Energy_ETD1_2D(N,xa,xb,tau,T,epsilon,u0)


Nx = N; Ny = N; nu = 1; Nt = round(T/tau); hx = (xb-xa)/Nx;  hy = hx; 
F = @(u) 1/4*(1-u.^2).^2; f = @(u) u.^3 - u;
beta = sqrt(4/3*(1+epsilon^2/nu)); kappa = 3;





kx = [0:Nx/2-1, -Nx/2:-1]; ky = [0:Ny/2-1, -Ny/2:-1];
[KX, KY] = meshgrid(kx, ky);

L_hat = -4/hx^2 * sin(pi*KX/Nx).^2 - 4/hy^2 * sin(pi*KY/Ny).^2; 
Lk_hat = epsilon^2/nu * L_hat - kappa;

z = tau*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) ( exp(z)-1 ) ./ z;
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u = u0; t = (0:Nt)'*tau; E = zeros(Nt+1,1); 
E(1) = Energy_FD_2D(u0,epsilon,F,hx,hy);
for n = 1:Nt
    u_hat = fft2(u); 
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0); 
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifft2(Nu_hat));

    Nu = max(min(Nu, kappa*beta), -kappa*beta); Nu_hat = fft2(Nu);
    
    u_new_hat = phi0(z) .* u_hat + tau * phi1(z) .* Nu_hat;
    u_new = real(ifft2(u_new_hat)); 
    u = u_new;

    E(n+1,1) = Energy_FD_2D(u,epsilon,F,hx,hy);
end
end


function E_FD = Energy_FD_2D(u,epsilon,F,hx,hy)


ux = (circshift(u, [0, -1]) - circshift(u, [0, 1])) / (2*hx);
uy = (circshift(u, [-1, 0]) - circshift(u, [1, 0])) / (2*hy);


E_FD = hx * hy * sum(sum( 1/2*epsilon^2*(ux.^2 + uy.^2) + F(u) ));
end