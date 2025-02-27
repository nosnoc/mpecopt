clear; clc; close all

%% Very good example where global optimum leads to more iterations, and reduced lpec is useful!
% Example  MPCC * requires a large penalty parametres
import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
a = 4;
f_sym = @(x,y) a*(x-1).^2+(y-1).^2;
f_grad_sym = @(x,y) [8*(x-1);2*(y-1)];

f = f_sym(x1,x2);
g = [];
x0 = [1;0];
% x0 = [0;1];

lbx = [0;0];
ubx = [inf;inf];
G = x1;
H = x2;
lbg = [];
ubg = [];
mpec = struct('x', x,'f',f, 'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
%%  Settings
solver_settings = MPECOptimizerOptions();
solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
solver_settings.consider_all_complementarities_in_lpec = 1;
solver_settings.rho_TR_phase_ii_init = 1.5;
solver_settings.TR_reducing_factor = 0.5;
solver_settings.compute_tnlp_stationary_point  = 0;
solver_settings.plot_lpec_iterate = 1;
solver_settings.stop_if_S_stationary = 0;

solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
x_opt = full(result_active_set.x);
f_opt = full(result_active_set.f);
fprintf('solution is (%2.4f,%2.4f) \n',x_opt(1),x_opt(2));
fprintf('objective is %2.4f \n',f_opt );
fprintf('\n');
%%
nice_plot_colors
tt = 0:1:5;
figure
plot(tt,tt*0,'k','LineWidth',1.5);
hold on
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
axis equal
plot(stats_active_set.iter.X_outer(1,:),stats_active_set.iter.X_outer(2,:),'rs')
plot(stats_active_set.iter.X_lpec(1,:),stats_active_set.iter.X_lpec(2,:),'bd')


if 1
    for i = 1:size(stats_active_set.iter.X_outer,2)-1
        quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
        text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
    end
else
    for i = 1:size(stats_active_set.iter.X_outer,2)-1
        quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), -2*(a*stats_active_set.iter.X_outer(1,i)-1),-2*(stats_active_set.iter.X_outer(2,i)-1), 0.5, 'Color',matlab_red,'LineWidth',2);
        text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
    end
end
xlim([-1 5])
ylim([-1 5])

% Generate grid points for x and y
x = linspace(-1, 5, 50);
y = linspace(-1, 5, 50);
% Create a grid of points
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f_sym(X, Y);
% Plot the contour lines
contour(X, Y, Z, 30);
xlabel('x');
ylabel('y');
%%
% x_opt = [1;0];
x_opt = [0;1];

figure
latexify_plot()
filename = 'global_local_mpec1.pdf';
linewidht = 2;
fontsize = 16;
markersize = 10;
rho = 1.2;
rho = 0.5;
% Given function
% Define grid
[x1, x2] = meshgrid(-1.2:0.01:3, -1.2:0.01:3);
nice_plot_colors
z = f_sym(x1, x2);
% Plot contour lines
contour(x1, x2, z,15,'LineWidth',linewidht-0.5);
hold on;
dx1 = 2*(a*x1-1);
dx2 = 2*(x2-1);
tt = 0:0.5:3;
f_grad = -f_grad_sym(x_opt(1),x_opt(2));
quiver(x_opt(1),x_opt(2),f_grad(1),f_grad(2),'LineWidth',linewidht+0.5,'AutoScaleFactor',0.2,'Color',matlab_red)

% case 1
quiver(x_opt(1),x_opt(2),0,0,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_magenta)
% case 2
% quiver(x_opt(1),x_opt(2),1.2,-1,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_magenta)
% case 3
% quiver(x_opt(1),x_opt(2),-1,1.2,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_magenta)


plot(tt*0,tt,'k','LineWidth',linewidht)
plot(tt,tt*0,'k','LineWidth',linewidht)
plot(1,0,'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
plot(0,1,'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)

% Trust region radius
x = x_opt;
% trust region radius
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
hold on
plot(tt,tt*0-rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');

x_bar = [1,0];
f_bar = ['$f(\bar{x}) = ' sprintf('%2.2f$',f_sym(x_bar(1),x_bar(2)))];
text(0.9,+0.3, '$\bar{x}$','FontSize',fontsize)
% text(1,-0.5, f_bar,'FontSize',fontsize)

x_hat = [0,1];
f_hat = ['$f(\hat{x}) = ' sprintf('%2.2f$',f_sym(x_hat(1),x_hat(2)))];
text(-0.25,1, '$\hat{x}$','FontSize',fontsize)
% text(-0.5,1.2, f_hat,'FontSize',fontsize)
% Latexify labels
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
axis equal
xlim([-1.2 3])
ylim([-1.2 3])
% legend({'$f(x)$', '$\nabla f(x^*)$', '$d$','$0\leq x_1 \perp x_2\geq 0$' },'BackgroundAlpha',0.9)
legend({'$f(x)$', '$-\nabla f(\bar{x})$', '$d$' },'BackgroundAlpha',0.9)
set(gca,'FontSize',fontsize)
exportgraphics(gcf, filename, 'ContentType', 'vector')







