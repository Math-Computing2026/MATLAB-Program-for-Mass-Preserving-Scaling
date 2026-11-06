
clear; clc; close all;
xa = 0; xb = 2*pi; N = 128; Nx = N; Ny = N; T = 100; tau = 1.5; epsilon = 0.01; 

u0 = 0.9575 * (2*rand(Ny,Nx)-1);
[Mass_Err1] = Mass_ETD1(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err2] = Mass_ETD2(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err3] = Mass_rETDRK3(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err4] = Mass_rETDRK4(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err5] = Mass_rETDRK5(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err6] = Mass_rETDRK6(N,xa,xb,tau,T,epsilon,u0);
Nt = round(T/tau); t = (1:Nt)' * tau;  





C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6]; markers = {'o','+','v','s','p','d'};
figure(1)  
semilogy(t,Mass_Err1,'--o','Color',C1,'LineWidth',1.2,'MarkerSize',6); hold on
semilogy(t,Mass_Err2,':x','Color',C2,'LineWidth',1.2,'MarkerSize',6); hold on
semilogy(t,Mass_Err3,'-v' ,'Color',C3,'LineWidth',1.2,'MarkerSize',6); hold on
semilogy(t,Mass_Err4,'--d', 'Color',C4,'LineWidth',1.2,'MarkerSize',6); hold on
semilogy(t,Mass_Err5,'-.p' ,'Color',C5,'LineWidth',1.2,'MarkerSize',6); hold on
semilogy(t,Mass_Err6,'-s' ,'Color',C6,'LineWidth',1.2,'MarkerSize',6); hold on
legend('ETD1','ETDRK2','rETDRK3','rETDRK4','rETDRK5','rETDRK6','FontSize', 12,'location','southeast');
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('Time','FontSize',16,'FontName','Times New Roman'); 
ylabel('Mass Error','FontSize',16,'FontName','Times New Roman');
xlim([0, T]); 
