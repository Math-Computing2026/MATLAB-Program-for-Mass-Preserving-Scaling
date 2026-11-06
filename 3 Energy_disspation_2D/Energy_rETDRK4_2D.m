

function [u,E,t] = Energy_rETDRK4_2D(N,xa,xb,tau,T,epsilon,u0)

Nx = N; Ny = N; nu = 1; Nt = round(T/tau); hx = (xb-xa)/Nx;  hy = hx; 
F = @(u) 1/4*(1-u.^2).^2; f = @(u) u.^3 - u; 
kappa = 3; beta = sqrt(4/3*(1+epsilon^2/nu));





kx = [0:Nx/2-1, -Nx/2:-1]; ky = [0:Ny/2-1, -Ny/2:-1];
[KX, KY] = meshgrid(kx, ky);

L_hat = -4/hx^2 * sin(pi*KX/Nx).^2 - 4/hy^2 * sin(pi*KY/Ny).^2;  
Lk_hat = epsilon^2/nu * L_hat - kappa;

z = tau*Lk_hat; z_half = tau/2*Lk_hat; z_13 = tau/3*Lk_hat; z_23 = 2*tau/3*Lk_hat; 
phi0 = @(z) exp(z);
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + 1/2.*(abs(z)<=1e-14);
phi3 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2)./(z.^3)) + 1/6.*(abs(z)<=1e-14);
phi4 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6)./(z.^4)) + 1/24.*(abs(z)<=1e-14);
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u = u0; t = (0:Nt)'*tau; E = zeros(Nt+1,1); 
E(1) = Energy_FD_2D(u0,epsilon,F,hx,hy);
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; 
    F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifft2(Nu_hat));
    

    w1_hat = phi0(z).*u_hat + tau*phi1(z).*Nu_hat;
    w1_half_hat = phi0(z_half).*u_hat + tau/2*phi1(z_half).*Nu_hat;
    w1_tau_3_hat = phi0(z_13) .* u_hat + tau/3*phi1(z_13) .* Nu_hat;
    w1_2tau_3_hat = phi0(z_23) .* u_hat + 2*tau/3*phi1(z_23) .* Nu_hat;
    w1 = real(ifft2(w1_hat)); 
    
    F1 = f(w1) - epsilon^2/nu*w1; F1_hat = fft2(F1);
    Nw1_hat = B_hat .* F1_hat + kappa * w1_hat;
 

    c1_1_hat = Nw1_hat - Nu_hat;  
    w2_hat = w1_hat + tau*phi2(z) .* c1_1_hat; 
    w2_half_hat = w1_half_hat + tau*(1/2)^2*phi2(z_half).*c1_1_hat;
    w2 = real(ifft2( w2_hat)); 
    w2_half = real(ifft2( w2_half_hat));
    
    F2_half = f(w2_half) - epsilon^2/nu*w2_half; 
    F2_half_hat = fft2(F2_half);
    F2 = f(w2) - epsilon^2/nu*w2; 
    F2_hat = fft2(F2);
    
    Nw2_s1_hat = B_hat .* F2_half_hat + kappa * w2_half_hat;
    Nw2_s2_hat = B_hat .* F2_hat + kappa * w2_hat;
    

    b1_hat = Nw2_s1_hat - Nu_hat; b2_hat = Nw2_s2_hat - Nu_hat;
    c2_1_hat = 4*b1_hat - b2_hat; c2_2_hat = -4*b1_hat + 2*b2_hat;
    
    w3_s1_hat = w1_tau_3_hat + tau*(1/3)^2*phi2(z_13).*c2_1_hat + 2*tau*(1/3)^3*phi3(z_13).*c2_2_hat;
    w3_s2_hat = w1_2tau_3_hat + tau*(2/3)^2*phi2(z_23).*c2_1_hat + 2*tau*(2/3)^3*phi3(z_23).*c2_2_hat;
    w3_s3_hat = w1_hat + tau*phi2(z).*c2_1_hat + 2*tau*phi3(z).*c2_2_hat;
    
    w3_s1 = real(ifft2(w3_s1_hat));  w3_s2 = real(ifft2(w3_s2_hat));  w3_s3 = real(ifft2(w3_s3_hat)); 

    F_13 = f(w3_s1) - epsilon^2/nu*w3_s1;  F_13_hat = fft2(F_13); 
    F_23 = f(w3_s2) - epsilon^2/nu*w3_s2;   F_23_hat = fft2(F_23); 
    F_3 = f(w3_s3) - epsilon^2/nu*w3_s3;   F_3_hat = fft2(F_3);
    
    Nw3_s1_hat = B_hat .* F_13_hat + kappa*w3_s1_hat;  
    Nw3_s2_hat = B_hat .* F_23_hat + kappa*w3_s2_hat;
    Nw3_s3_hat = B_hat .* F_3_hat + kappa*w3_s3_hat;  
     
    d1_hat = Nw3_s1_hat - Nu_hat; d2_hat = Nw3_s2_hat - Nu_hat; d3_hat = Nw3_s3_hat - Nu_hat;
    
    c3_1_hat = 9*d1_hat - 9/2*d2_hat + d3_hat;     c3_1 = real(ifft2(c3_1_hat));         
    c3_2_hat = -45/2*d1_hat + 18*d2_hat - 9/2*d3_hat;    c3_2 = real(ifft2(c3_2_hat));  
    c3_3_hat = 27/2*d1_hat - 27/2*d2_hat + 9/2*d3_hat;   c3_3 = real(ifft2(c3_3_hat));    



    xi = linspace(0,1,51); xi = xi(2:end);
    alpha = 1;
    for k = 1:length(xi)
        R_s = xi(k)*c3_1 + xi(k)^2*c3_2 + xi(k)^3*c3_3;
        idx = find(abs(R_s)>1e-14);
        if ~isempty(idx)
            local_alpha = (kappa*beta - sign(R_s(idx)) .* Nu(idx)) ./ abs(R_s(idx));
            alpha = min(alpha, min(local_alpha));
        end
    end
    alpha = min(alpha,1); alpha = max(alpha,0);
    u_new_hat = phi0(z) .* u_hat + tau*phi1(z) .* Nu_hat...
                + alpha *( tau*phi2(z).*c3_1_hat + 2*tau*phi3(z).*c3_2_hat + 6*tau*phi4(z).*c3_3_hat);
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