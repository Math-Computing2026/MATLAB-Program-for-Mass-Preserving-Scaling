
clear all; clc; close all;
N = 128; xa = -0.4; xb = 0.4; epsilon = 0.01; nu = 1;
beta = sqrt(4/3*(1+epsilon^2/nu)); tau = 3; tau_ref = 0.1; T = 400;  
hx = (xb - xa)/N; x = xa + (0:N-1)*hx; y = x; 
[X, Y] = meshgrid(x, y); 


X1 = -0.15; Y1 = 0.0; R1 = 0.16; X2 = 0.15;  Y2 = 0.0; R2 = 0.12;

r1 = sqrt((X - X1).^2 + (Y - Y1).^2);
r2 = sqrt((X - X2).^2 + (Y - Y2).^2);
u0 = 1 - tanh((r1 - R1)/(sqrt(2)*epsilon)) - tanh((r2 - R2)/(sqrt(2)*epsilon));
[u2,E2,t] = Energy_ETDRK2_2D(N,xa,xb,tau,T,epsilon,u0); 
[u3,E3,t] = Energy_rETDRK3_2D(N,xa,xb,tau,T,epsilon,u0); 
[u5,E5,t] = Energy_rETDRK5_2D(N,xa,xb,tau,T,epsilon,u0); 
[u6,E6,t] = Energy_rETDRK6_2D(N,xa,xb,tau,T,epsilon,u0); 
[u_ref,E_ref,t_ref] = Energy_rETDRK4_2D(N,xa,xb,tau_ref,T,epsilon,u0); 
Nt = round(T/tau); Nt_ref = round(T/tau_ref);


C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6]; markers = {'o','+','v','s','p','d'};
figure(1)
hold on; box on;
labels = {'ETDRK2','rETDRK3','rETDRK5','rETDRK6','Reference'};
plot(t,E2,'-.o','Color',C1,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:4:Nt);  hold on 
plot(t,E3,'--p','Color',C5,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:4:Nt);  
plot(t,E5,'-v','Color',C3,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:4:Nt);  
plot(t,E6,'-d','Color',C4,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:4:Nt); 
plot(t_ref,E_ref,':x', 'Color',C2,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:55:Nt_ref); hold on
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('Time','FontSize',16,'FontName','Times New Roman');
ylabel('Energy','FontSize',16,'FontName','Times New Roman');
legend(labels, 'Location', 'northeast', 'Fontsize',12);
xlim([0,T]);

