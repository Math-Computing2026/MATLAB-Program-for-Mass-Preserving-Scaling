
clear all; clc; close all;
N = 128; xa = -0.80; xb = 0.80; nu = 1; epsilon = 0.02;  T = 120;      
hx = (xb-xa)/N; x = xa + (0:N-1)*hx; y = x; z = x; tau = 1; 
[X, Y, Z] = ndgrid(x, y, z); 
beta = sqrt(4/3*(1+epsilon^2/nu));


u0 = 0.2 * tanh(100*(0.2 - abs(X-0.15).^0.8 - abs(Y).^0.8 - abs(Z-0.15).^0.8))...
     + 0.2 * tanh(100*(0.2 - abs(X+0.15).^0.8 - abs(Y).^0.8 - abs(Z-0.15).^0.8))...
     + 0.2 * tanh(100*(0.2 - abs(X).^0.8 - abs(Y).^0.8 - abs(Z+0.1).^0.8)) + 0.3;
[u2,u2_hist,t] = MBP_ETDRK2_3D(N,xa,xb,tau,T,epsilon,u0);
[u6,u6_hist,t] = MBP_rETDRK6_3D(N,xa,xb,tau,T,epsilon,u0);
Nt = round(T/tau); 

C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];

line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6];
markers = {'o','.','v','s','p','d'};
figure(1)
plot(t,u2_hist,'--o' ,'Color',C1,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:2:Nt);  hold on
plot(t,u6_hist,'-d' ,'Color',C4,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:2:Nt);  hold on
set(gca,'FontSize',12,'FontName','Times New Roman'); box on;
yline(beta, 'k--', 'LineWidth', 1.8, 'Alpha', 0.8);
xlabel('Time','FontSize',16,'FontName','Times New Roman'); xlim([0, T]); 
ylabel('Supermum norm','FontSize',16,'FontName','Times New Roman');
legend('ETDRK2','rETDRK6', '$\beta$', 'Interpreter', 'latex','Location', 'Southeast', 'FontSize', 12);




inset_ax = axes('Position', [0.53, 0.47, 0.24, 0.22]); 
hold on; box on;
plot(t, u2_hist, '--o',  'Color', C1, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:3:Nt); 
plot(t, u6_hist, '-d',  'Color', C4, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:3:Nt); 

 xlim(inset_ax, [100, 120]); ylim(inset_ax, [0.99, 1.02]); 
set(inset_ax, 'FontSize', 10, 'FontName', 'Times New Roman');