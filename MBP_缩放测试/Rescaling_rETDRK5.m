

function [u,u_hist,t] = Rescaling_rETDRK5(N,xa,xb,tau,T,epsilon,u0)

nu = 1; Nt = round(T/tau); Nx = N; Ny = N;
hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; hy = hx; y = x;
[X,Y] = meshgrid(x,y);

theta = 0.8; theta_c = 1.6; u_safe = @(u) max(min(u, 1-1e-6), -1+1e-6);  
f = @(u) theta/2*log((1+u_safe(u))./(1-u_safe(u))) - theta_c*u; 
beta = 0.957504024077268;  kappa = 10; 

kx = [0:Nx/2-1,-Nx/2:-1]'; ky = [0:Ny/2-1,-Ny/2:-1]';
Lx = -4/hx^2 * sin(pi*kx/Nx).^2; Ly = -4/hy^2 * sin(pi*ky/Ny).^2; 
[LX,LY] = meshgrid(Lx,Ly);
L_hat =  LX + LY;
Lk_hat = epsilon^2/nu * L_hat - kappa;


z = tau*Lk_hat; z_half = tau/2*Lk_hat; z_13 = tau/3*Lk_hat; z_23 = 2*tau/3*Lk_hat; 
z_14 = tau/4*Lk_hat; z_34 = 3*tau/4*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + 1/2.*(abs(z)<=1e-14);
phi3 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2)./(z.^3)) + 1/6.*(abs(z)<=1e-14);
phi4 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6)./(z.^4)) + 1/24.*(abs(z)<=1e-14);
phi5 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6-z.^4/24)./(z.^5)) + 1/120.*(abs(z)<=1e-14);
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);


u = u0; t = (0:Nt)'*tau;
u_hist = zeros(Nt+1,1); u_hist(1) = max(abs(u0(:)));
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat; 
    Nu = real(ifft2(Nu_hat));
    




    Nu = max(min(Nu, kappa*beta), -kappa*beta);
    Nu_hat = fft2(Nu);



    w1_hat = phi0(z) .* u_hat + tau*phi1(z) .* Nu_hat;                  
    w1_half_hat = phi0(z_half).* u_hat + tau/2*phi1(z_half).* Nu_hat;    
    w1_tau_3_hat = phi0(z_13) .* u_hat + tau/3*phi1(z_13) .* Nu_hat;    
    w1_2tau_3_hat = phi0(z_23) .* u_hat + 2*tau/3*phi1(z_23) .* Nu_hat;    
    w1_tau_4_hat = phi0(z_14) .* u_hat + tau/4*phi1(z_14) .* Nu_hat;      
    w1_3tau_4_hat = phi0(z_34) .* u_hat + 3*tau/4*phi1(z_34) .* Nu_hat;   
    
    w1 = real(ifft2(w1_hat)); 
    
    F1 = f(w1) - epsilon^2/nu*w1; F1_hat = fft2(F1);
    Nw1_hat = B_hat .* F1_hat + kappa * w1_hat;
 

    c1_1_hat = Nw1_hat - Nu_hat;  
    w2_hat = w1_hat + tau*phi2(z) .* c1_1_hat; 
    w2_half_hat = w1_half_hat + tau*(1/2)^2*phi2(z_half).*c1_1_hat;
    w2 = real(ifft2( w2_hat)); w2_half = real(ifft2( w2_half_hat));
    
    F2_half = f(w2_half) - epsilon^2/nu*w2_half; F2_half_hat = fft2(F2_half);
    F2 = f(w2) - epsilon^2/nu*w2; F2_hat = fft2(F2);
    
    Nw2_s1_hat = B_hat .* F2_half_hat + kappa * w2_half_hat;
    Nw2_s2_hat = B_hat .* F2_hat + kappa * w2_hat;
    

    b1_hat = Nw2_s1_hat - Nu_hat; b2_hat = Nw2_s2_hat - Nu_hat;
    c2_1_hat = 4*b1_hat - b2_hat; c2_2_hat = -4*b1_hat + 2*b2_hat;
    
    w3_s1_hat = w1_tau_3_hat + tau*(1/3)^2*phi2(z_13).*c2_1_hat + 2*tau*(1/3)^3*phi3(z_13).*c2_2_hat;
    w3_s2_hat = w1_2tau_3_hat + tau*(2/3)^2*phi2(z_23).*c2_1_hat + 2*tau*(2/3)^3*phi3(z_23).*c2_2_hat;
    w3_s3_hat = w1_hat + tau*phi2(z).*c2_1_hat + 2*tau*phi3(z).*c2_2_hat;
    w3_s1 = real(ifft2(w3_s1_hat)); w3_s2 = real(ifft2(w3_s2_hat)); w3_s3 = real(ifft2(w3_s3_hat)); 
    

    F_13 = f(w3_s1) - epsilon^2/nu*w3_s1; 
    F_23 = f(w3_s2) - epsilon^2/nu*w3_s2; 
    F_3 = f(w3_s3) - epsilon^2/nu*w3_s3;
    F_13_hat = fft2(F_13); F_23_hat = fft2(F_23); F_3_hat = fft2(F_3);
    Nw3_s1_hat = B_hat .* F_13_hat + kappa*w3_s1_hat;  
    Nw3_s2_hat = B_hat .* F_23_hat + kappa*w3_s2_hat;
    Nw3_s3_hat = B_hat .* F_3_hat + kappa*w3_s3_hat;  
     
    d1_hat = Nw3_s1_hat - Nu_hat; d2_hat = Nw3_s2_hat - Nu_hat; d3_hat = Nw3_s3_hat - Nu_hat;
    
    c3_1_hat = 9*d1_hat - 9/2*d2_hat + d3_hat;                
    c3_2_hat = -45/2*d1_hat + 18*d2_hat - 9/2*d3_hat;       
    c3_3_hat = 27/2*d1_hat - 27/2*d2_hat + 9/2*d3_hat;       
    
    w4_s1_hat = w1_tau_4_hat + tau*(1/4)^2*phi2(z_14).*c3_1_hat...
                + 2*tau*(1/4)^3*phi3(z_14).*c3_2_hat + 6*tau*(1/4)^4*phi4(z_14).*c3_3_hat;
    w4_s2_hat = w1_half_hat + tau*(1/2)^2*phi2(z_half).*c3_1_hat...
                + 2*tau*(1/2)^3*phi3(z_half).*c3_2_hat + 6*tau*(1/2)^4*phi4(z_half).*c3_3_hat;
    w4_s3_hat = w1_3tau_4_hat + tau*(3/4)^2*phi2(z_34).*c3_1_hat...
                + 2*tau*(3/4)^3*phi3(z_34).*c3_2_hat + 6*tau*(3/4)^4*phi4(z_34).*c3_3_hat;
    w4_s4_hat = w1_hat + tau*phi2(z).*c3_1_hat + 2*tau*phi3(z).*c3_2_hat + 6*tau*phi4(z).*c3_3_hat;
    w4_s1 = real(ifft2(w4_s1_hat));  w4_s2 = real(ifft2(w4_s2_hat));
    w4_s3 = real(ifft2(w4_s3_hat));  w4_s4 = real(ifft2(w4_s4_hat));
    

    F_14 = f(w4_s1) - epsilon^2/nu*w4_s1;  F_14_hat = fft2(F_14);
    F_24 = f(w4_s2) - epsilon^2/nu*w4_s2;  F_24_hat = fft2(F_24);
    F_34 = f(w4_s3) - epsilon^2/nu*w4_s3;  F_34_hat = fft2(F_34);
    F_4 = f(w4_s4) - epsilon^2/nu*w4_s4;   F_4_hat = fft2(F_4);
    
    Nw4_s1_hat =  B_hat .* F_14_hat + kappa*w4_s1_hat;  
    Nw4_s2_hat =  B_hat .* F_24_hat + kappa*w4_s2_hat;  
    Nw4_s3_hat =  B_hat .* F_34_hat + kappa*w4_s3_hat;  
    Nw4_s4_hat =  B_hat .* F_4_hat + kappa*w4_s4_hat;  
    
    e1_hat = Nw4_s1_hat - Nu_hat; e2_hat = Nw4_s2_hat - Nu_hat;
    e3_hat = Nw4_s3_hat - Nu_hat; e4_hat = Nw4_s4_hat - Nu_hat;
    
    c4_1_hat = 16*e1_hat - 12*e2_hat + 16/3*e3_hat - e4_hat;  c4_1 = real(ifft2(c4_1_hat));
    c4_2_hat = -208/3*e1_hat + 76*e2_hat - 112/3*e3_hat + 22/3*e4_hat; c4_2 = real(ifft2(c4_2_hat));
    c4_3_hat = 96*e1_hat - 128*e2_hat + 224/3*e3_hat - 16*e4_hat; c4_3 = real(ifft2(c4_3_hat));
    c4_4_hat = -128/3*e1_hat + 64*e2_hat - 128/3*e3_hat + 32/3*e4_hat; c4_4 = real(ifft2(c4_4_hat));
    

    xi = linspace(0,1,21); xi = xi(2:end);
    alpha = 1;
    for k = 1:length(xi)
        R_s = xi(k)*c4_1 + xi(k)^2*c4_2 + xi(k)^3*c4_3 + xi(k)^4*c4_4;  
        idx = find(abs(R_s)>1e-14);
        if ~isempty(idx)
            local_alpha = (kappa*beta - sign(R_s(idx)) .* Nu(idx)) ./ abs(R_s(idx));
            alpha = min(alpha, min(local_alpha));
        end
    end
    alpha = min(alpha,1); alpha = max(alpha,0); 
    
    u_new_hat = w1_hat + alpha * (tau*phi2(z).*c4_1_hat + 2*tau*phi3(z).*c4_2_hat...
                + 6*tau*phi4(z).*c4_3_hat + 24*tau*phi5(z).*c4_4_hat);
    u_new = real(ifft2(u_new_hat));
    u = u_new; u_hist(n+1,1) = max(abs(u(:)));
end
end