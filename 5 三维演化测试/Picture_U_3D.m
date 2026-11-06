
clear; clc; close all;
N = 128; xa = -0.80; xb = 0.80; nu = 1; epsilon = 0.02;  T = 2;      
hx = (xb-xa)/N; x = xa + (0:N-1)*hx; y = x; z = x; tau = 0.2; 
[X, Y, Z] = meshgrid(x, y, z);
u0 = 0.2 * tanh(100*(0.2 - abs(X-0.15).^0.8 - abs(Y).^0.8 - abs(Z-0.15).^0.8))...
     + 0.2 * tanh(100*(0.2 - abs(X+0.15).^0.8 - abs(Y).^0.8 - abs(Z-0.15).^0.8))...
     + 0.2 * tanh(100*(0.2 - abs(X).^0.8 - abs(Y).^0.8 - abs(Z+0.1).^0.8)) + 0.3;

[u2] = U_ETDRK2_3D(N,xa,xb,tau,T,epsilon,u0,nu);
[u6] = U_rETDRK6_3D(N,xa,xb,tau,T,epsilon,u0,nu);

C1 = [0.00,0.80,0.75]; C2 = [0.45, 0.25, 0.60]; 
figure(1)
fv1 = isosurface(X, Y, Z, u2, 0.5); 
p1 = patch(fv1); p1.FaceColor = C1; p1.EdgeColor = 'none';
p1.FaceAlpha = 0.8; p1.BackFaceLighting = 'lit';
daspect([1 1 1]); view(3); view(-30,25);       
axis tight; camlight; lighting gouraud;  
xlim([xa,xb]); ylim([xa,xb]); zlim([xa,xb]);
set(gca,'FontSize',12,'FontName','Times New Roman','LineWidth', 1); 
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);
figure(2)
fv2 = isosurface(X, Y, Z, u6, 0.5); 
p2 = patch(fv2); p2.FaceColor = C1; p2.EdgeColor = 'none'; 
p2.FaceAlpha = 0.8; p2.BackFaceLighting = 'lit';
daspect([1 1 1]); view(3); view(-30,25);          
axis tight; camlight; lighting gouraud; 
xlim([xa,xb]); ylim([xa,xb]); zlim([xa,xb]); 
set(gca,'FontSize',12,'FontName','Times New Roman','LineWidth', 1);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);



