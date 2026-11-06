
function [u,u_hist,t] = Rescaling_rETDRK3(N,xa,xb,tau,T,epsilon,u0)


nu = 1; Nt = round(T/tau); Nx = N; Ny = N;
hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; hy = hx; y = x;
[X,Y] = meshgrid(x,y);
f = @(u) u.^3 - u; beta = 1; kappa = 3;




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


u = u0; t = (0:Nt)'*tau;
u_hist = zeros(Nt+1,1); u_hist(1) = max(abs(u0(:)));
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;  Nu = real(ifft2(Nu_hat)); 
    Nu = max(min(Nu, kappa*beta), -kappa*beta); 
    Nu_hat = fft2(Nu);  

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


    xi = linspace(0,1,51); xi = xi(2:end);
    alpha = 1;
    for k = 1:length(xi)
        R_s = xi(k)*c2_1 + xi(k)^2*c2_2;
        idx = find(abs(R_s)>1e-14);
        if ~isempty(idx)
            local_alpha = (kappa*beta - sign(R_s(idx)) .* Nu(idx)) ./ abs(R_s(idx));
            alpha = min(alpha, min(local_alpha));
        end
    end
    alpha = min(alpha,1); alpha = max(alpha,0);
    
    u_new_hat = phi0(z) .* u_hat + tau*phi1(z) .* Nu_hat...
                + alpha * (tau*phi2(z).*c2_1_hat + 2*tau*phi3(z).*c2_2_hat);
    u_new = real(ifft2(u_new_hat));
    u = u_new; u_hist(n+1,1) = max(abs(u(:)));
end
end

