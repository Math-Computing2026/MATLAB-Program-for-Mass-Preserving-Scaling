


















clear all; clc; close all; epsilon = 0.01;
xa = 0; xb = 2*pi; tau = 2; nu = 1; 
N = 128; Nx = N; Ny = N; T = 100;

beta = 0.957504024077268; u0 = 0.9575 * (2*rand(Ny,Nx)-1);
[u1,u1_hist,t] = Standard_ETDRK5(N,xa,xb,tau,T,epsilon,u0);
[u2,u2_hist,t] = Rescaling_rETDRK5(N,xa,xb,tau,T,epsilon,u0);

C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C6 = [0.15,0.15,0.15];
figure; 
plot(t,u1_hist,':','Color',C2,'Marker','none','LineWidth',2.3); hold on;
plot(t,u2_hist,'-.','Color',C1,'Marker','none','LineWidth',2);  hold on;
yline(beta,'--','Color',C6,'LineWidth',2,'Alpha',0.8); box on;

set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
xlabel('Time', 'FontSize', 16,'FontName', 'Times New Roman');
ylabel('Supremum norm', 'FontSize', 16,'FontName', 'Times New Roman');
legend('ETDRK5','rETDRK5','$\beta$','Interpreter','latex','Location','southeast');



inset_ax = axes('Position', [0.30, 0.25, 0.24, 0.22]); hold on;
plot(t,u1_hist,':','Color',C2,'Marker','none','LineWidth',2.3);
plot(t,u2_hist,'-.','Color',C1,'Marker','none','LineWidth',2);
yline(beta,'--','Color',C6,'LineWidth',2,'Alpha',0.8); box on;

xlim(inset_ax, [90, 100]); ylim(inset_ax, [0.957501, 0.957505]);
set(inset_ax, 'FontSize', 10, 'FontName', 'Times New Roman');







