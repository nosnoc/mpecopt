%%
clear all
clc
close all
import casadi.*

%% Example  MPCC - Is unbounded for small penalties
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
p = SX.sym('p');
f = x1^2+x2^2-4*x1*x2;
x0 = [0.5;0.6];
x0 = [0.1;4];
x0 = [3;1];
G = x1;
H = x2;
lbx = [0;0];
ubx = [inf;inf];

g = [];
lbg = [];
ubg = [];
mpec = struct('x', x, 'f', f, 'g', g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% Scholtes
settings = HomotopySolverOptions();
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
x_opt_homotopy = full(result_homotopy.x);
fprintf('x_opt = (%2.4f,%2.4f), f_opt = %2.4f. \n',x_opt_homotopy(1),x_opt_homotopy(2),f_opt_homotopy);
%%  Settings
solver_settings = mpecopt.Options();
solver_settings.settings_lpec.lpec_solver ="Gurobi";
% solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
% solver_settings.initialization_strategy = "RelaxAndProject";
solver_settings.consider_all_complementarities_in_lpec = true;
solver_settings.tol_B_stationarity = 1e-8;
solver_settings.relax_and_project_iters = 2;
solver_settings.relax_and_project_kappa = 0.5;
solver_settings.relax_and_project_sigma0 = 0.01;
solver_settings.plot_lpec_iterate = 1;
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% [result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);


solver = mpecopt.Solver(mpec, solver_settings);
[result_active_set,stats_active_set] = solver.solve(solver_initalization);

x_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);


fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Scholtes \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('\n');
fprintf(' || x_reg - x_active_set || = %2.2e \n',norm(x_opt_homotopy-x_opt_active_set));

fprintf('solution is (%2.4f,%2.4f) \n',x_opt_active_set(1),x_opt_active_set(2));

%%
nice_plot_colors
figure
xline(0,'k');
hold on
grid on
yline(0,'k');
axis equal
plot(stats_active_set.iter.X_outer(1,:),stats_active_set.iter.X_outer(2,:),'rs')
for i = 1:size(stats_active_set.iter.X_outer,2)-1
    quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
end
f = @(x,y) x.^2+y.^2-4*x*y;
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
