
function [Mass_Err] = Mass_rETDRK6(N,xa,xb,tau,T,epsilon,u0)


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

z = tau*Lk_hat; z_half = tau/2*Lk_hat; z_13 = tau/3*Lk_hat; z_23 = 2*tau/3*Lk_hat; 
z_14 = tau/4*Lk_hat; z_34 = 3*tau/4*Lk_hat; 
z_15 = tau/5*Lk_hat; z_25 = 2*tau/5*Lk_hat; z_35 = 3*tau/5*Lk_hat; z_45 = 4*tau/5*Lk_hat;
phi0 = @(z) exp(z);
phi1 = @(z) (abs(z)>1e-14).*((exp(z)-1)./z) + (abs(z)<=1e-14).*1;
phi2 = @(z) (abs(z)>1e-14).*((exp(z)-1-z)./(z.^2)) + 1/2.*(abs(z)<=1e-14);
phi3 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2)./(z.^3)) + 1/6.*(abs(z)<=1e-14);
phi4 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6)./(z.^4)) + 1/24.*(abs(z)<=1e-14);
phi5 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6-z.^4/24)./(z.^5)) + 1/120.*(abs(z)<=1e-14);
phi6 = @(z) (abs(z)>1e-14).*((exp(z)-1-z-z.^2/2-z.^3/6-z.^4/24-z.^5/120)./(z.^6)) + 1/720.*(abs(z)<=1e-14);
B_hat = 1/nu * (1 ./ (1-nu*L_hat) - 1);

u = u0; t = (0:Nt)'*tau;
t = (0:Nt)'*tau; Mass = zeros(Nt+1,1); Mass_Err = zeros(Nt,1);  
Mass(1) = hx*hy*sum(u0(:));
for n = 1:Nt
    u_hat = fft2(u);
    
    F0 = f(u) - epsilon^2/nu*u; F0_hat = fft2(F0);
    Nu_hat = B_hat .* F0_hat + kappa * u_hat;
    Nu = real(ifft2(Nu_hat));

    w1_hat = phi0(z) .* u_hat + tau*phi1(z) .* Nu_hat;
    w1_half_hat = phi0(z_half).* u_hat + tau/2*phi1(z_half).* Nu_hat;
    w1_tau_3_hat = phi0(z_13) .* u_hat + tau/3*phi1(z_13) .* Nu_hat;
    w1_2tau_3_hat = phi0(z_23) .* u_hat + 2*tau/3*phi1(z_23) .* Nu_hat;
    w1_tau_4_hat = phi0(z_14) .* u_hat + tau/4*phi1(z_14) .* Nu_hat;
    w1_3tau_4_hat = phi0(z_34) .* u_hat + 3*tau/4*phi1(z_34) .* Nu_hat;
    w1_tau_5_hat = phi0(z_15) .* u_hat + tau/5*phi1(z_15) .* Nu_hat;
    w1_2tau_5_hat = phi0(z_25) .* u_hat + 2*tau/5*phi1(z_25) .* Nu_hat;
    w1_3tau_5_hat = phi0(z_35) .* u_hat + 3*tau/5*phi1(z_35) .* Nu_hat;
    w1_4tau_5_hat = phi0(z_45) .* u_hat + 4*tau/5*phi1(z_45) .* Nu_hat;
    
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
    
    c4_1_hat = 16*e1_hat - 12*e2_hat + 16/3*e3_hat - e4_hat;
    c4_2_hat = -208/3*e1_hat + 76*e2_hat - 112/3*e3_hat + 22/3*e4_hat;
    c4_3_hat = 96*e1_hat - 128*e2_hat + 224/3*e3_hat - 16*e4_hat;
    c4_4_hat = -128/3*e1_hat + 64*e2_hat - 128/3*e3_hat + 32/3*e4_hat;
    
    w5_s1_hat = w1_tau_5_hat + tau*(1/5)^2*phi2(z_15).*c4_1_hat + 2*tau*(1/5)^3*phi3(z_15).*c4_2_hat...
            + 6*tau*(1/5)^4*phi4(z_15).*c4_3_hat + 24*tau*(1/5)^5*phi5(z_15).*c4_4_hat;
    w5_s2_hat = w1_2tau_5_hat + tau*(2/5)^2*phi2(z_25).*c4_1_hat + 2*tau*(2/5)^3*phi3(z_25).*c4_2_hat...
            + 6*tau*(2/5)^4*phi4(z_25).*c4_3_hat + 24*tau*(2/5)^5*phi5(z_25).*c4_4_hat;
    w5_s3_hat = w1_3tau_5_hat + tau*(3/5)^2*phi2(z_35).*c4_1_hat + 2*tau*(3/5)^3*phi3(z_35).*c4_2_hat...
            + 6*tau*(3/5)^4*phi4(z_35).*c4_3_hat + 24*tau*(3/5)^5*phi5(z_35).*c4_4_hat; 
    w5_s4_hat = w1_4tau_5_hat + tau*(4/5)^2*phi2(z_45).*c4_1_hat + 2*tau*(4/5)^3*phi3(z_45).*c4_2_hat...
            + 6*tau*(4/5)^4*phi4(z_45).*c4_3_hat + 24*tau*(4/5)^5*phi5(z_45).*c4_4_hat;
    w5_s5_hat = w1_hat + tau*phi2(z).*c4_1_hat + 2*tau*phi3(z).*c4_2_hat...
                   + 6*tau*phi4(z).*c4_3_hat + 24*tau*phi5(z).*c4_4_hat;
    w5_s1 = real(ifft2(w5_s1_hat)); w5_s2 = real(ifft2(w5_s2_hat)); w5_s3 = real(ifft2(w5_s3_hat));
    w5_s4 = real(ifft2(w5_s4_hat)); w5_s5 = real(ifft2(w5_s5_hat));

    F_15 = f(w5_s1) - epsilon^2/nu*w5_s1; F_15_hat = fft2(F_15);
    F_25 = f(w5_s2) - epsilon^2/nu*w5_s2; F_25_hat = fft2(F_25);
    F_35 = f(w5_s3) - epsilon^2/nu*w5_s3; F_35_hat = fft2(F_35);
    F_45 = f(w5_s4) - epsilon^2/nu*w5_s4; F_45_hat = fft2(F_45);
    F_5 = f(w5_s5) - epsilon^2/nu*w5_s5;  F_5_hat = fft2(F_5);
    
    Nw5_s1_hat =  B_hat .* F_15_hat + kappa*w5_s1_hat;  
    Nw5_s2_hat =  B_hat .* F_25_hat + kappa*w5_s2_hat;  
    Nw5_s3_hat =  B_hat .* F_35_hat + kappa*w5_s3_hat;   
    Nw5_s4_hat =  B_hat .* F_45_hat + kappa*w5_s4_hat;   
    Nw5_s5_hat =  B_hat .* F_5_hat + kappa*w5_s5_hat;  
    
    g1_hat = Nw5_s1_hat - Nu_hat; g2_hat = Nw5_s2_hat - Nu_hat;
    g3_hat = Nw5_s3_hat - Nu_hat; g4_hat = Nw5_s4_hat - Nu_hat; g5_hat = Nw5_s5_hat - Nu_hat;
    
    c5_1_hat = 25*g1_hat - 25*g2_hat + 50/3*g3_hat - 25/4*g4_hat + g5_hat;
    c5_2_hat = -1925/12*g1_hat + 2675/12*g2_hat - 325/2*g3_hat + 1525/24*g4_hat - 125/12*g5_hat;
    c5_3_hat = 8875/24*g1_hat - 7375/12*g2_hat + 6125/12*g3_hat - 5125/24*g4_hat + 875/24*g5_hat;
    c5_4_hat = -4375/12*g1_hat + 8125/12*g2_hat - 625*g3_hat + 6875/24*g4_hat - 625/12*g5_hat;
    c5_5_hat = 3125/24*g1_hat - 3125/12*g2_hat + 3125/12*g3_hat - 3125/24*g4_hat + 625/24*g5_hat;
    c5_1 = real(ifft2(c5_1_hat)); c5_2 = real(ifft2(c5_2_hat)); c5_3 = real(ifft2(c5_3_hat)); 
    c5_4 = real(ifft2(c5_4_hat)); c5_5 = real(ifft2(c5_5_hat)); 


    xi = linspace(0,1,51); xi = xi(2:end);
    alpha = 1;
    for k = 1:length(xi)
        R_s = xi(k)*c5_1 + xi(k)^2*c5_2 + xi(k)^3*c5_3 + xi(k)^4*c5_4 + xi(k)^5*c5_5;
        idx = find(abs(R_s)>1e-14);
        if ~isempty(idx)
            local_alpha = (kappa*beta - sign(R_s(idx)) .* Nu(idx)) ./ abs(R_s(idx));
            alpha = min(alpha, min(local_alpha));
        end
    end
    alpha = min(alpha,1); alpha = max(alpha,0);
 
    u_new_hat = w1_hat + alpha * ( tau*phi2(z).*c5_1_hat + 2*tau*phi3(z).*c5_2_hat...
                + 6*tau*phi4(z).*c5_3_hat + 24*tau*phi5(z).*c5_4_hat + 120*tau*phi6(z).*c5_5_hat );
    u_new = real(ifft2(u_new_hat));
    u = u_new;
    Mass(n+1,1) = hx*hy*sum(u(:)); Mass_Err(n,1) = abs(Mass(n+1) - Mass(1));
end
end