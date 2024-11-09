clear; clc; close all

%% Projection lpec fail unless sigma is quite small
import casadi.*
a = 100;
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
a = 1.1;
tau0 = 1;
f_sym = @(x,y) -2*x+1*y;
g_sym  = @(x,y) -x.^2-(y-a).^2+1; % good lpec inconsistent bnlp
g_nabla_sym = @(x,y)[-2*x;-2*(y-a)];

% g_sym  = @(x,y) y-(x).^2; % hard to verify B
f = f_sym(x1,x2);
g = [g_sym(x1,x2)];
x0 = [0.0;0.5];

lbx = [0;0];
ubx = [inf;1];
G = x1;
H = x2;
lbg = 0;
ubg = inf;

mpec = struct('x', x,'f',f, 'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% homotopy
settings = HomotopySolverOptions();
settings.tol  = 1e-12;
settings.max_iter = 3;
settings.kappa = 0.2;
settings.homotopy_parameter_steering = 'Direct';
settings.sigma0 = tau0;
settings.check_B_stationarity  = 0;
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
x_opt_homotopy = full(result_homotopy.x);
fprintf('x_homotopy = (%2.4f,%2.4f), f_opt = %2.4f. \n',x_opt_homotopy(1),x_opt_homotopy(2),f_opt_homotopy);

tau0 = (settings.kappa)^(settings.max_iter-1)*tau0;
%%
solver_settings = MPECOptimizerOptions();
% solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% [result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
% x_opt_active_set = full(result_active_set.x);
% f_opt_active_set = full(result_active_set.f);
% fprintf('x_active_set = (%2.4f,%2.4f), f_opt = %2.4f.  \n',x_opt_active_set(1),x_opt_active_set(2),f_opt_active_set);
%%
% x_star = [1;1];
filename1  = 'phase1_mpec_2.pdf';
filename2  = 'phase1_lpec_2.pdf';
nice_plot_colors;
linewidht = 2;
fontsize = 16;
markersize = 10;
x_opt = x_opt_homotopy;
nice_plot_colors
tt = 0:1:5;
figure
latexify_plot();
x = linspace(-2, 5, 50);
y = linspace(-2, 5, 50);
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f_sym(X, Y);
% Plot the contour lines
contour(X, Y, Z, 30);
hold on;
fimplicit(g_sym, [-3 3], 'MeshDensity', 100, 'LineWidth',2,'Color',matlab_blue)
fimplicit(@(x,y) x*y-tau0, [0 10], 'MeshDensity', 100, 'LineWidth',2,'Color',matlab_green)
plot(x_opt(1),x_opt(2),'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
text(x_opt(1)-0.2,x_opt(2)+0.3, '${x}^*(\tau)$','FontSize',fontsize)

plot(tt,tt*0,'k','LineWidth',1.5);
hold on
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
axis equal
xlim([-1.2 4])
ylim([-1.2 4])
xlabel('$x_1$');
ylabel('$x_2$');

set(gca,'FontSize',fontsize);
legend({'$f(x)$', '$c(x)=0$', '$x_1x_2=\tau$' },'BackgroundAlpha',0.9)
exportgraphics(gcf, filename1, 'ContentType', 'vector')

% Linearization
figure
latexify_plot();
contour(X, Y, Z, 30);
hold on
quiver(x_opt(1),x_opt(2),2,-1,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_red)
fimplicit(g_sym, [-3 3], 'MeshDensity', 100, 'LineWidth',1, 'LineStyle','--','Color',matlab_blue)
fimplicit(@(x,y) g_sym(x_opt(1),x_opt(2))+g_nabla_sym(x_opt(1),x_opt(2))'*([x;y]-x_opt), [-5 5], 'MeshDensity', 100, 'LineWidth',linewidht,'Color',matlab_blue)
% optimal lpec solution
% quiver(x_opt(1),x_opt(2),-x_opt(1)+0.5,-x_opt(2)-0.1,'LineWidth',linewidht+0.5,'Color',matlab_magenta,'MaxHeadSize',0.8)
quiver(x_opt(1),x_opt(2),-x_opt(1),-x_opt(2)+0.08,'LineWidth',linewidht+0.5,'Color',matlab_magenta,'MaxHeadSize',0.8)



plot(x_opt(1),x_opt(2),'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
text(x_opt(1)-0.2,x_opt(2)+0.3, '${x}^*(\tau)$','FontSize',fontsize)


% Trust region radius
rho = 1;
x = x_opt;
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
hold on
plot(tt,tt*0-rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');


tt = 0:1:5;
plot(tt,tt*0,'k','LineWidth',1.5);
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
axis equal
xlim([-1.2 4])
ylim([-1.2 4])
legend({'$f(x)$', '$-\nabla f(x^*(\tau))$', '$c(x)=0$', '$c(x^*(\tau))\!+\! \nabla c(x^*(\tau))^\top d=0$', '$d$' },'BackgroundAlpha',0.8,'NumColumns',1,'FontSize',fontsize-2)
xlabel('$x_1$');
ylabel('$x_2$');
set(gca,'FontSize',fontsize);

exportgraphics(gcf, filename2, 'ContentType', 'vector')

%%
% figure
