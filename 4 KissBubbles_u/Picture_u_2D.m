
clear; clc; close all;
N = 128; xa = -0.4; xb = 0.4; epsilon = 0.01; nu = 0.8;
tau = 3; tau_ref = 0.01; T = 500;  
hx = (xb - xa)/N; x = xa + (0:N-1)*hx; y = x; 
[X, Y] = meshgrid(x, y); 


X1 = -0.15; Y1 = 0.0; R1 = 0.16; X2 = 0.15;  Y2 = 0.0; R2 = 0.12;

r1 = sqrt((X - X1).^2 + (Y - Y1).^2);
r2 = sqrt((X - X2).^2 + (Y - Y2).^2);
u0 = 1 - tanh((r1 - R1)/(sqrt(2)*epsilon)) - tanh((r2 - R2)/(sqrt(2)*epsilon));
u0 = max(min(u0, 1), -1);


[U1,t] = U_ETDRK2_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U3,t] = U_rETDRK3_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U5,t] = U_rETDRK5_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U6,t] = U_rETDRK6_2D(N,xa,xb,tau,T,epsilon,u0,nu);
[U_ref,t] = U_rETDRK4_2D(N,xa,xb,tau_ref,T,epsilon,u0,nu);

x = xa + (0:N)*hx; y = x; [X, Y] = meshgrid(x, y); 
u1 = zeros(N+1,N+1); u1(1:N,1:N) = U1; u1(N+1,1:N) = U1(1,:);u1(1:N,N+1) = U1(:,1);  u1(end,end) = U1(1,1);
u3 = zeros(N+1,N+1); u3(1:N,1:N) = U3; u3(N+1,1:N) = U3(1,:); u3(1:N,N+1) = U3(:,1); u3(end,end) = U3(1,1);  
u5 = zeros(N+1,N+1); u5(1:N,1:N) = U5; u5(N+1,1:N) = U5(1,:); u5(1:N,N+1) = U5(:,1); u5(end,end) = U5(1,1);
u6 = zeros(N+1,N+1); u6(1:N,1:N) = U6; u6(N+1,1:N) = U6(1,:);  u6(1:N,N+1) = U6(:,1); u6(end,end) = U6(1,1); 
u_ref = zeros(N+1,N+1); u_ref(1:N,1:N) = U_ref; u_ref(N+1,1:N) = U_ref(1,:); 
u_ref(1:N,N+1) = U_ref(:,1); u_ref(end,end) = U_ref(1,1);

figure(1)
contourf(X, Y, u1, 60, 'LineColor', 'none'); hold on;
axis square; colormap('jet');  colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
axis([-0.4,0.4,-0.4,0.4]); xticks(-0.4:0.2:0.4);  yticks(-0.4:0.2:0.4);  

figure(2)
contourf(X, Y, u3, 60, 'LineColor', 'none'); hold on;
axis square; colormap('jet');  colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
axis([-0.4,0.4,-0.4,0.4]); xticks(-0.4:0.2:0.4);  yticks(-0.4:0.2:0.4);  

figure(3)
contourf(X, Y, u5, 60, 'LineColor', 'none'); hold on;
axis square; colormap('jet');  colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
axis([-0.4,0.4,-0.4,0.4]); xticks(-0.4:0.2:0.4);  yticks(-0.4:0.2:0.4);  

figure(4)
contourf(X, Y, u6, 60, 'LineColor', 'none'); hold on;
axis square; colormap('jet');  colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
axis([-0.4,0.4,-0.4,0.4]); xticks(-0.4:0.2:0.4);  yticks(-0.4:0.2:0.4);  

figure(5)
contourf(X, Y, u_ref, 60, 'LineColor', 'none'); hold on;
axis square; colormap('jet'); colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
axis([-0.4,0.4,-0.4,0.4]); xticks(-0.4:0.2:0.4);  yticks(-0.4:0.2:0.4);  