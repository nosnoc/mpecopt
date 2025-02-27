clear; clc; close all

%% Very good example where global optimum leads to more iterations, and reduced lpec is useful!
% Example  MPCC * requires a large penalty parametres
import casadi.*
a = 10;
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
a = 100;
f = (a*x1-1)^2+(x2-1)^2;
g = [];
x0 = [0.0104;1.25];
x0 = [2;0];
x0 = [0;4];

lbx = [0;0];
ubx = [inf;inf];
G = x1;
H = x2;
lbg = [];
ubg = [];

mpec = struct('x', x,'f',f, 'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% homotopy
settings = HomotopySolverOptions();
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);
fprintf('x_opt = (%2.4f,%2.4f), f_opt = %2.4f. \n',w_opt_homotopy(1),w_opt_homotopy(2),f_opt_homotopy);
%%  Settings
% x0 = [1/a;0];
x0 = [1/a;0];
solver_settings = MPECOptimizerOptions();
solver_settings.settings_lpec.lpec_solver ="Gurobi";
solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
solver_settings.consider_all_complementarities_in_lpec = 1;
solver_settings.rho_TR_phase_ii_init = 1;
solver_settings.TR_reducing_factor = 0.5;
solver_settings.tol_B_stationarity = 1e-6;
solver_settings.plot_mpec_multipliers = 0;
solver_settings.plot_lpec_iterate = 1;


solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);
fprintf('solution is (%2.4f,%2.4f) \n',w_opt_active_set(1),w_opt_active_set(2));


fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('homotopy \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
% fprintf('Objectives: homotopy: %2.2e, Active set: %2.2e \n',f_opt,f_opt_active_set)
% fprintf('Biactive: homotopy: %d, Active set: %d \n',stats_homotopy.n_biacive,stats_active_set.n_biacive)
fprintf('\n');
fprintf(' || w_homotopy - w_active_set || = %2.2e \n',norm(w_opt_homotopy-w_opt_active_set));

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

f = @(x,y) (a*x-1).^2+(y-1).^2;
% Generate grid points for x and y
x = linspace(-1, 5, 50);
y = linspace(-1, 5, 50);
% Create a grid of points
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f(X, Y);
% Plot the contour lines
contour(X, Y, Z, 30);
xlabel('x');
ylabel('y');
%%

figure
latexify_plot()
a = 4;
filename = 'lpec_small_rho1.pdf';
linewidht = 2;
fontsize = 16;
markersize = 10;
rho = 1.1;
% Given function
f = @(x1, x2) (a*x1 - 1).^2 + (x2 - 1).^2;
% Define grid
[x1, x2] = meshgrid(-1.2:0.01:2, -1.2:0.01:2);

% Evaluate function at grid points
nice_plot_colors
z = f(x1, x2);
% Plot contour lines
contour(x1, x2, z,15,'LineWidth',linewidht-0.5);
hold on;
dx1 = 2*(a*x1-1);
dx2 = 2*(x2-1);
tt = 0:0.5:2;
quiver(1/a,0,0,1,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_red)

if rho < 0.9/a
    quiver(1/a,0,0,0,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_magenta)
else
    quiver(1/a,0,-1/a,rho,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_magenta)
end

plot(tt*0,tt,'k','LineWidth',linewidht)
plot(tt,tt*0,'k','LineWidth',linewidht)
plot(1/a,0,'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
plot(0,1,'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)

% Trust region radius
x = [1/a,0];
% trust region radius
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
hold on
plot(tt,tt*0-rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
text(1/a,-0.2, '$\bar{x}$','FontSize',fontsize)
text(-0.2,1, '$\hat{x}$','FontSize',fontsize)
% Latexify labels
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
axis equal
xlim([-1.2 2])
ylim([-1.2 2])
% legend({'$f(x)$', '$\nabla f(x^*)$', '$d$','$0\leq x_1 \perp x_2\geq 0$' },'BackgroundAlpha',0.9)
legend({'$f(x)$', '$-\nabla f(\bar{x})$', '$d$' },'BackgroundAlpha',0.9)
set(gca,'FontSize',fontsize)
exportgraphics(gcf, filename, 'ContentType', 'vector')







