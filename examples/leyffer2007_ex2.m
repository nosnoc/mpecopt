clear all; clc; close all;
import casadi.*
%% Leyffer 2007, example 2, origin is C stationary, (1,0) is S-stationary
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
% parameters
p = SX.sym('p');
p0 = 0;

f = (x1-1)^2+x2^3+x2^2+p;
G = x1;
H = x2;

% x0 = [3;0];
x0 = [0;2];
% x0 = [0;0];
lbx = [0;0];
ubx = [inf;inf];
g = [2*x1+x2];
lbg = [0.1];
ubg = [inf];

mpec = struct('x', x, 'f', f, 'g', g,'G',G,'H',H,'p',p);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg,'p0',p0);
% Scholtes
settings = HomotopySolverOptions();
[result_scholtes,stats_scholtes] = mpec_homotopy_solver(mpec,solver_initalization,settings);
w_opt = full(result_scholtes.x);
f_opt_scholtes = full(result_scholtes.f);
%%
% Pivoting
solver_settings = MPECOptimizerOptions();
solver_settings.settings_lpec.lpec_solver ="Highs";
solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
solver_settings.consider_all_complementarities_in_lpec = 1; % indentifying the biactive set fails and wrong sol?
% solver_settings.tol_active = 1e-5;
% solver_settings.consider_all_complementarities_in_lpec = 1; % 
solver_settings.plot_lpec_iterate = 1;
solver_settings.tol_B_stationarity = 1e-8;
solver_settings.rho_TR_phase_ii_init = 0.3;
solver_settings.rho_TR_phase_i_init = 10;
solver_settings.stop_if_S_stationary = 0;

[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);

fprintf('solution is (%2.4f,%2.4f) \n',w_opt_active_set(1),w_opt_active_set(2));

% solution.X_outer

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
    text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
end

f = @(x,y) (x-1).^2 + y.^2 + y.^3;
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
