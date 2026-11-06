
function [u,t] = U_ETDRK2_3D(N,xa,xb,tau,T,epsilon,u0,nu)

Nx = N; Ny = N; Nz = N;
Nt = round(T/tau); hx = (xb-xa)/Nx; hy = hx; hz = hx;
f = @(u) u.^3-u; beta = sqrt(4/3*(1+epsilon^2/nu)); kappa = 3;





kx = [0:Nx/2-1,-Nx/2:-1]; ky = [0:Ny/2-1,-Ny/2:-1]; kz = [0:Nz/2-1,-Nz/2:-1];
[Kx,Ky,Kz] = ndgrid(kx,ky,kz); 
LX = -4/hx^2 * sin(pi*Kx/Nx).^2; LY = -4/hy^2 * sin(pi*Ky/Ny).^2; LZ = -4/hz^2 * sin(pi*Kz/Nz).^2; 
L_hat = LX + LY + LZ;
Lk_hat = epsilon^2/nu * L_hat - kappa;

z_tau = tau*Lk_hat;
phi0 = @(z) exp(z);
phi1=@(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14);
phi2=@(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + (abs(z)<=1e-14).*0.5;

B_hat = 1/nu * (1./(1 - nu*L_hat)-1);

u = u0; t=(0:Nt)'*tau;
for n = 1:Nt
    u_hat = fftn(u);

    F0 = f(u) - epsilon^2/nu*u; F0_hat = fftn(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifftn(Nu_hat));

    w1_hat = phi0(z_tau).*u_hat  + tau*phi1(z_tau).*Nu_hat;
    w1 = real(ifftn(w1_hat));
    N1 = f(w1) - epsilon^2/nu*w1; N1_hat = fftn(N1);
    Nw1_hat = B_hat.*N1_hat + kappa*w1_hat;


    c1_1_hat = Nw1_hat - Nu_hat;
    u_new_hat = w1_hat + tau*phi2(z_tau).*c1_1_hat;
    u_new = real(ifftn(u_new_hat)); u = u_new;
end

end