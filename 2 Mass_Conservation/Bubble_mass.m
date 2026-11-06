
clear all; clc; close all;
N = 256; xa = -0.4; xb = 0.4; epsilon = 0.01; nu = 1;
beta = sqrt(4/3*(1+epsilon^2/nu)); tau = 3; tau_ref = 0.1; T = 500;  
hx = (xb - xa)/N; x = xa + (0:N-1)*hx; y = x; 
[X, Y] = meshgrid(x, y); 


X1 = -0.15; Y1 = 0.0; R1 = 0.16; X2 = 0.15;  Y2 = 0.0; R2 = 0.12;

r1 = sqrt((X - X1).^2 + (Y - Y1).^2);
r2 = sqrt((X - X2).^2 + (Y - Y2).^2);
u0 = 1 - tanh((r1 - R1)/(sqrt(2)*epsilon)) - tanh((r2 - R2)/(sqrt(2)*epsilon));
 
[Mass_Err2] = Mass_ETD2(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err3] = Mass_rETDRK3(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err5] = Mass_rETDRK5(N,xa,xb,tau,T,epsilon,u0);
[Mass_Err6] = Mass_rETDRK6(N,xa,xb,tau,T,epsilon,u0);
[MassErr_ref] = Mass_rETDRK4(N,xa,xb,tau_ref,T,epsilon,u0);
Nt = round(T/tau); t = (1:3:Nt)' * tau;  Nt_ref = round(T/tau_ref); t_ref = (1:45:Nt_ref)'*tau_ref;
Mass_Err2 = Mass_Err2(1:3:Nt); Mass_Err3 = Mass_Err3(1:3:Nt); Mass_Err5 = Mass_Err5(1:3:Nt); 
Mass_Err6 = Mass_Err6(1:3:Nt); MassErr_ref = MassErr_ref(1:45:Nt_ref);

C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6]; markers = {'o','+','v','s','p','d'};

figure(1)  
semilogy(t,Mass_Err2,'-.o','Color',C1,'LineWidth',1.3,'MarkerSize',6 ); hold on
semilogy(t,Mass_Err3,'--p' ,'Color',C5,'LineWidth',1.3,'MarkerSize',6 ); 
semilogy(t,Mass_Err5,'-.v' ,'Color',C3,'LineWidth',1.3,'MarkerSize',6 ); 
semilogy(t,Mass_Err6,'-d' ,'Color',C4,'LineWidth',1.3,'MarkerSize',6 ); 
semilogy(t_ref,MassErr_ref,':x', 'Color',C2,'LineWidth',1.3,'MarkerSize',6 ); 
legend('ETDRK2','rETDRK3','rETDRK5','rETDRK6','Reference','FontSize', 12,'location','southeast');
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('Time','FontSize',16,'FontName','Times New Roman'); 
ylabel('Mass Error','FontSize',16,'FontName','Times New Roman');
xlim([0, T]); 
