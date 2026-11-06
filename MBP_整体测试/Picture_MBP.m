
clear all; clc; close all;
xa = 0; xb = 2*pi; tau = 1.5; nu = 1; epsilon = 0.01;
N = 128; Nx = N; Ny = N; T = 100; 
beta = sqrt(4/3*(1+epsilon^2/nu)); u0 = 0.01*rand(Ny,Nx)-0.05;

[u1,u1_hist,t] = Standard_ETD1(N,xa,xb,tau,T,epsilon,u0);
[u2,u2_hist,t] = Standard_ETDRK2(N,xa,xb,tau,T,epsilon,u0);
[u3,u3_hist,t] = Rescaling_rETDRK3(N,xa,xb,tau,T,epsilon,u0);
[u4,u4_hist,t] = Rescaling_rETDRK4(N,xa,xb,tau,T,epsilon,u0);
[u5,u5_hist,t] = Rescaling_rETDRK5(N,xa,xb,tau,T,epsilon,u0);
[u6,u6_hist,t] = Rescaling_rETDRK6(N,xa,xb,tau,T,epsilon,u0);




C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
figure(1)

line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6];
markers = {'o','.','v','s','p','d'};
u_hist = {u1_hist, u2_hist, u3_hist, u4_hist, u5_hist, u6_hist};
labels = {'ETD1', 'ETDRK2', 'rETDRK3', 'rETDRK4', 'rETDRK5', 'rETDRK6'};

plot(t,u1_hist,'--o','Color',C1,'LineWidth',1.3,'MarkerSize', 6); hold on
plot(t,u2_hist,':x' ,'Color',C2,'LineWidth',1.3,'MarkerSize', 6);  hold on 
plot(t,u3_hist,'-v' ,'Color',C3,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,u4_hist,'--d', 'Color',C4,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,u5_hist,'-.p' ,'Color',C5,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,u6_hist,'-s' ,'Color',C6,'LineWidth',1.3,'MarkerSize', 6);  hold on
set(gca,'FontSize',12,'FontName','Times New Roman'); box on;
yline(beta, 'k--', 'LineWidth', 1.8, 'Alpha', 0.8);
xlabel('Time','FontSize',16,'FontName','Times New Roman');
ylabel('Supermum norm','FontSize',16,'FontName','Times New Roman');
legend([labels, '$\beta$'], 'Interpreter', 'latex','Location', 'Southeast', 'FontSize', 12);
xlim([0, T]); 