function [] = plot_lpec(nabla_f, x, d, rho_TR)
nice_plot_colors
figure
grid on
axis equal
% trust region radius
tt = linspace(x(1)-rho_TR,x(1)+rho_TR,3);
plot(tt,tt*0+rho_TR+x(2),'color',matlab_blue,'LineWidth',2,'HandleVisibility','off');
hold on
plot(tt,tt*0-rho_TR+x(2),'color',matlab_blue,'LineWidth',2,'HandleVisibility','off');
tt = linspace(x(2)-rho_TR,x(2)+rho_TR,3);
plot(tt*0+x(1)+rho_TR,tt,'color',matlab_blue,'LineWidth',2,'HandleVisibility','off');
plot(tt*0+x(1)-rho_TR,tt,'color',matlab_blue,'LineWidth',2,'HandleVisibility','off');
% comp set 
tt = 0:0.1:max(5,max(x(:))+rho_TR+1);
plot(tt,tt*0,'k','LineWidth',2,'HandleVisibility','off');
hold on
plot(tt*0,tt,'k','LineWidth',2,'HandleVisibility','off');

% important points and vectors
x_trail = x + d;
plot(x(1),x(2),'r.','MarkerSize',25,'LineWidth',2,'DisplayName','$x^k$');
hold on
plot(x_trail (1),x_trail (2),'kh','MarkerSize',5,'LineWidth',2,'DisplayName','$x^k+d$');
quiver(x(1), x(2), d(1),d(2), 1, 'Color',matlab_red,'LineWidth',1.0,'DisplayName','$d$');
quiver(x(1), x(2), -nabla_f(1)/norm(nabla_f), -nabla_f(2)/norm(nabla_f), 1, 'Color',matlab_magenta,'LineWidth',1.0,'DisplayName','$-\nabla f(x^k)$');
legend('Interpreter','latex')


% level curves
f = @(x,y) nabla_f(1)*x+nabla_f(2)*y;
X = linspace(-1, 5, 50);
Y = linspace(-1, 5, 50);
[X, Y] = meshgrid(X, Y);
Z = f(X, Y);
contour(X, Y, Z, 30,'HandleVisibility','off');

axis equal
% xline(x(1)+rho_TR,'color',matlab_blue,'LineWidth',2)
% xline(x(1)-rho_TR,'color',matlab_blue,'LineWidth',2)
% yline(x(2)+rho_TR,'color',matlab_blue,'LineWidth',2)
% yline(x(2)-rho_TR,'color',matlab_blue,'LineWidth',2)

end