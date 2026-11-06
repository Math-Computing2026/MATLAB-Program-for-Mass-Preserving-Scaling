
clear all; clc; close all;
Nx = 512; Ny = Nx; T = 1; kappa = 3; hx = 2*pi/Nx; hy = 2*pi/Ny;
[U] = Refer_ETDRK6(Nx,Ny,0.1*2^(-8),T,kappa);
Tau = 2.^(-(1:8));
for i = 1:length(Tau)
   tau = Tau(i);








   [u3] = Rescaling_ETDRK3_VCH(Nx,Ny,tau,T,kappa);
   MaxErr3(i) = max(max(abs(U-u3)));
   L2Err3(i) = sqrt(hx*hy*sum(sum((U-u3).^2)));




   



     



end




C3 = [0.00,0.00,1.00];
C4 = [1.00,0.00,0.00]; C5 = [1.00,0.07,0.65];
color_Ref = [0.25, 0.25, 0.25];



























figure(3);
loglog(Tau,L2Err3,'-o','Color',C4,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor','w'); hold on;
loglog(Tau,MaxErr3,'-^','Color',C3,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor','w');

C3 = L2Err3(1) / (Tau(1)^3) * 15; 
loglog(Tau,C3*Tau.^3,'--','Color',color_Ref,'LineWidth',2); grid on; 
set(gca,'GridLineStyle',':','GridAlpha',0.5); box on;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('Time step size','FontSize',16,'FontName','Times New Roman');
ylabel('Error','FontSize',16,'FontName','Times New Roman');
lgd = legend('$L^2$ norm error','$L^\infty$ norm error','3th-order','Interpreter','latex','Location','southeast');
set(lgd,'FontSize',14,'Box','on');








































