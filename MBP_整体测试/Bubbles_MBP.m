
clear all; clc; close all;
N = 128; xa = -0.4; xb = 0.4; epsilon = 0.01; nu = 1;
beta = sqrt(4/3*(1+epsilon^2/nu)); tau = 3; tau_ref = 0.1; T = 400;  
hx = (xb - xa)/N; x = xa + (0:N-1)*hx; y = x; 
[X, Y] = meshgrid(x, y); 


X1 = -0.15; Y1 = 0.0; R1 = 0.16; X2 = 0.15;  Y2 = 0.0; R2 = 0.12;

r1 = sqrt((X - X1).^2 + (Y - Y1).^2);
r2 = sqrt((X - X2).^2 + (Y - Y2).^2);
u0 = 1 - tanh((r1 - R1)/(sqrt(2)*epsilon)) - tanh((r2 - R2)/(sqrt(2)*epsilon));

[u2,u2_hist,t] = Standard_ETDRK2(N,xa,xb,tau,T,epsilon,u0);
[u3,u3_hist,t] = Rescaling_rETDRK3(N,xa,xb,tau,T,epsilon,u0);
[u5,u5_hist,t] = Rescaling_rETDRK5(N,xa,xb,tau,T,epsilon,u0);
[u6,u6_hist,t] = Rescaling_rETDRK6(N,xa,xb,tau,T,epsilon,u0);
[u_ref,ref_hist,t_ref] = Rescaling_rETDRK4(N,xa,xb,tau_ref,T,epsilon,u0);
Nt = round(T/tau); Nt_ref = round(T/tau_ref);

C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];

line_styles = {'--',':', '-', '--', '-.', '-'}; 
colors = [C1; C2; C3; C4; C5; C6];
markers = {'o','.','v','s','p','d'};
u_hist = {u2_hist,u3_hist,u5_hist,u6_hist,ref_hist};
labels = {'ETDRK2','rETDRK3','rETDRK5','rETDRK6','Reference'};
figure(1)
plot(t,u2_hist,'-.o' ,'Color',C1,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:4:Nt);  hold on 
plot(t,u3_hist,'--p' ,'Color',C5,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:5:Nt);  hold on
plot(t,u5_hist,'-v' ,'Color',C3,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:6:Nt);  hold on
plot(t,u6_hist,'-d' ,'Color',C4,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:7:Nt);  hold on
plot(t_ref,ref_hist,':x', 'Color',C2,'LineWidth',1.3,'MarkerSize',6,'MarkerIndices',1:85:Nt_ref);  hold on
set(gca,'FontSize',12,'FontName','Times New Roman'); box on;
yline(beta, 'k--', 'LineWidth', 1.8, 'Alpha', 0.8);
xlabel('Time','FontSize',16,'FontName','Times New Roman');
ylabel('Supermum norm','FontSize',16,'FontName','Times New Roman');
legend([labels, '$\beta$'], 'Interpreter', 'latex','Location', 'best', 'FontSize', 12);
xlim([0, T]); 



inset_ax = axes('Position', [0.58, 0.30, 0.24, 0.22]); 
hold on; box on;
plot(t, u2_hist, '-.o', 'Color', C1, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:4:Nt); 
plot(t, u3_hist, '--p', 'Color', C5, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:5:Nt); 
plot(t, u5_hist, '-v',  'Color', C3, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:6:Nt); 
plot(t, u6_hist, '-d',  'Color', C4, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:7:Nt); 
plot(t_ref, ref_hist, ':x', 'Color', C2, 'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerIndices', 1:85:Nt_ref); 
xlim(inset_ax, [250, 320]); ylim(inset_ax, [0.99, 1.01]);  box on;
set(inset_ax, 'FontSize', 10, 'FontName', 'Times New Roman');