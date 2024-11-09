clear all
clc
close all
import casadi.*

%% Example MPCC
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];

f = x1+x2;
g = [x2^2-1; x1*x2];
x0 = [0.1;0.9];
G = x1;
H = x2;
lbx = [0;0];
ubx = [inf;inf];
lbg = [];
g = [];
ubg = [];

mpec = struct('x', x,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
settings = HomotopySolverOptions();
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);
%%  Settings
solver_settings = MPECOptimizerOptions();
solver_settings.consider_all_complementarities_in_lpec = false;
solver_settings.plot_lpec_iterate = 1;
% solver_settings.settings_lpec.lpec_solver = 'Ell_inf';

solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('homotopy \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf(' || w_homotopy - w_active_set || = %2.2e \n',norm(w_opt_homotopy-w_opt_active_set));
fprintf('solution is (%2.4f,%2.4f) \n',w_opt_active_set(1),w_opt_active_set(2));
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
for i = 1:size(stats_active_set.iter.X_outer,2)-1
    quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
end
f = @(x,y) x+y;
yline(1)
% Generate grid points for x and y
x = linspace(-1, 5, 50);
y = linspace(-1, 5, 50);
% Create a grid of points
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f(X, Y);
% Plot the contour lines
contour(X, Y, Z, 40);
xlabel('x');