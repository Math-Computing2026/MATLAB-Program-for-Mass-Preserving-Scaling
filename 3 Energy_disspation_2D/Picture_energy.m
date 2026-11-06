
clear all; clc; close all;
N = 128; Nx = N; Ny = N; tau = 1.5; nu = 1; epsilon = 0.01; T = 100;
xa = 0; xb = 2*pi; hx = (xb-xa)/Nx; j = (0:Nx-1)'; x = xa + j*hx; 
u0 = 0.01*rand(Ny,Nx)-0.05;



[E1,E1,t] = Energy_ETD1_2D(N,xa,xb,tau,T,epsilon,u0); 
[E2,E2,t] = Energy_ETDRK2_2D(N,xa,xb,tau,T,epsilon,u0); 
[u3,E3,t] = Energy_rETDRK3_2D(N,xa,xb,tau,T,epsilon,u0); 
[u4,E4,t] = Energy_rETDRK4_2D(N,xa,xb,tau,T,epsilon,u0); 
[u5,E5,t] = Energy_rETDRK5_2D(N,xa,xb,tau,T,epsilon,u0); 
[u6,E6,t] = Energy_rETDRK6_2D(N,xa,xb,tau,T,epsilon,u0); 




C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6]; markers = {'o','+','v','s','p','d'};
figure   
hold on; box on;
labels = {'ETD1', 'ETDRK2', 'rETDRK3', 'rETDRK4', 'rETDRK5', 'rETDRK6'};
plot(t,E1,'--o','Color',C1,'LineWidth',1.3,'MarkerSize', 6); hold on
plot(t,E2,':x' ,'Color',C2,'LineWidth',1.3,'MarkerSize', 6);  hold on 
plot(t,E3,'-v' ,'Color',C3,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,E4,'--d', 'Color',C4,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,E5,'-.p' ,'Color',C5,'LineWidth',1.3,'MarkerSize', 6);  hold on
plot(t,E6,'-s' ,'Color',C6,'LineWidth',1.3,'MarkerSize', 6);  hold on
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('Time','FontSize',16,'FontName','Times New Roman');
ylabel('Energy','FontSize',16,'FontName','Times New Roman');
legend(labels, 'Location', 'northeast', 'Fontsize',12);
xlim([0,T]);

