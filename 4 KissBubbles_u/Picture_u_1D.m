
clear; clc; close all;
N = 128; xa = -0.4; xb = 0.4; epsilon = 0.01; nu = 1; T = 420;        
hx = (xb - xa)/N; x = xa + (0:N-1)*hx; y = x; tau = 1; tau_ref = 0.01;
[X, Y] = meshgrid(x, y); 


X1 = -0.15; Y1 = 0.0; R1 = 0.16; 
X2 = 0.15;  Y2 = 0.0; R2 = 0.12;

r1 = sqrt((X - X1).^2 + (Y - Y1).^2);
r2 = sqrt((X - X2).^2 + (Y - Y2).^2);
u0 = 1 - tanh((r1 - R1)/(sqrt(2)*epsilon)) - tanh((r2 - R2)/(sqrt(2)*epsilon));
u0 = max(min(u0, 1), -1); 


[U2,t] = U_ETDRK2_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U3,t] = U_rETDRK3_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U5,t] = U_rETDRK5_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U6,t] = U_rETDRK6_2D(N,xa,xb,tau,T,epsilon,u0,nu); 
[U_ref,t] = U_rETDRK4_2D(N,xa,xb,tau_ref,T,epsilon,u0,nu);

x = xa + (0:N)*hx; y = x; [X, Y] = meshgrid(x, y); 
u2 = zeros(N+1,N+1); u2(1:N,1:N) = U2; u2(N+1,1:N) = U2(1,:);u2(1:N,N+1) = U2(:,1);  u2(end,end) = U2(1,1);
u3 = zeros(N+1,N+1); u3(1:N,1:N) = U3; u3(N+1,1:N) = U3(1,:); u3(1:N,N+1) = U3(:,1); u3(end,end) = U3(1,1);  
u5 = zeros(N+1,N+1); u5(1:N,1:N) = U5; u5(N+1,1:N) = U5(1,:); u5(1:N,N+1) = U5(:,1); u5(end,end) = U5(1,1);
u6 = zeros(N+1,N+1); u6(1:N,1:N) = U6; u6(N+1,1:N) = U6(1,:);  u6(1:N,N+1) = U6(:,1); u6(end,end) = U6(1,1); 
u_ref = zeros(N+1,N+1); u_ref(1:N,1:N) = U_ref; u_ref(N+1,1:N) = U_ref(1,:); 
u_ref(1:N,N+1) = U_ref(:,1); u_ref(end,end) = U_ref(1,1);

mid_idx = find(abs(y) == min(abs(y)), 1);  y_cut = y(mid_idx);
slice2 = u2(mid_idx,:);  slice3 = u3(mid_idx,:); slice5 = u5(mid_idx,:);
slice6 = u6(mid_idx,:); slice4 = u_ref(mid_idx,:);

C1 = [0.00, 0.45, 0.74];
C2 = [1.00,0.00,0.00];
C3 = [0.93,0.69,0.13];
C4 = [0.72, 0.27, 1.00];
C5 = [0.45, 0.74, 0.07];
C6 = [0.07, 0.62, 1.00];
figure(1)
plot(x, slice2, 'Color', C1, 'LineStyle', '-',  'LineWidth', 1.7); hold on;
plot(x, slice3, 'Color', C5, 'LineStyle', '-.', 'LineWidth', 1.7); hold on;
plot(x, slice5, 'Color', C3, 'LineStyle', '--', 'LineWidth', 1.7); hold on;
plot(x, slice6, 'Color', C4, 'LineStyle', '-.',  'LineWidth', 1.7);  hold on;
plot(x, slice4, 'Color', C2, 'LineStyle', ':',  'LineWidth', 1.8); hold on;
axis([xa xb -1.2 1.2]); 
set(gca,'FontSize',12,'FontName','Times New Roman'); box on;
xlabel('$x$','Interpreter','latex','FontSize',16); 
ylabel(['$u(x,0)$'],'Interpreter','latex','FontSize',16);
legend('ETDRK2','rETDRK3','rETDRK5','rETDRK6','Reference','Location','best','FontSize',12);


inset_ax = axes('Position', [0.38, 0.24, 0.24, 0.22]); hold on; 
plot(x,slice2,'-','Color',C1, 'LineWidth',1.7);
plot(x,slice3,'-.','Color',C5, 'LineWidth',1.7); 
plot(x,slice5,'-.','Color',C3, 'LineWidth',1.7); 
plot(x,slice6,'-.','Color',C4, 'LineWidth',1.7); 
plot(x,slice4,':','Color',C2, 'LineWidth',1.7); 
xlim(inset_ax, [0.192, 0.215]); ylim(inset_ax, [-0.98, -0.968]);  box on;
set(inset_ax, 'FontSize', 10, 'FontName', 'Times New Roman');
 
