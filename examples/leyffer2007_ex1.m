clear all; clc; close all;
import casadi.*
%% Leyffer 2007, example 2, origin is C stationary, (0,1) and (1,0) are S-stationary
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];

% parameters
p = SX.sym('p');
p0 = 2;

f = x1^2+x2^2-p*(x1+x2)+2;
G = x1;
H = x2;

x0 = [1;0.9];
% x0 = [1;0];
x0 = [0.5;1.5];
x0 = [0;0];
% x0 = [3;0];

lbx = [0;0];
ubx = [inf;inf];
g = []; lbg = []; ubg = [];

mpec = struct('x', x, 'f', f, 'g', g,'G',G,'H',H,'p',p);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg,'p0',p0);
% Scholtes
settings = HomotopySolverOptions();
[result_scholtes,stats_scholtes] = mpec_homotopy_solver(mpec,solver_initalization,settings);
w_opt = full(result_scholtes.x);
f_opt_scholtes = full(result_scholtes.f);
% Pivoting
solver_settings = mpecopt.Options();
solver_settings.settings_lpec.lpec_solver = "Highs_casadi";
solver_settings.initialization_strategy ="RelaxAndProject";
solver_settings.relax_and_project_homotopy_parameter_steering = "Ell_inf";
solver_settings.consider_all_complementarities_in_lpec = true;

% solver_settings.initialization_strategy ="TakeProvidedActiveSet";
% solver_initalization.y0 = 0;

solver = Mpecopt(mpec, solver_settings);
[result_active_set,stats_active_set] = solver.solve(solver_initalization);

w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);

fprintf('solution is (%2.4f,%2.4f) \n',w_opt_active_set(1),w_opt_active_set(2));
%%
nice_plot_colors
figure
title('Iterations of the active set method')
xline(0,'k');
hold on
grid on
yline(0,'k');
axis equal
plot(stats_active_set.iter.X_outer(1,:),stats_active_set.iter.X_outer(2,:),'rs')
for i = 1:size(stats_active_set.iter.X_outer,2)-1
    quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
    text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
end
f = @(x,y) x.^2+y.^2-2*(x+y)+2;
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
title('Iterations of the relaxation method')
xline(0,'k');
hold on
grid on
yline(0,'k');
axis equal
plot(stats_scholtes.iter.X_outer(1,:),stats_scholtes.iter.X_outer(2,:),'rs')
for i = 1:size(stats_scholtes.iter.X_outer,2)-1
    quiver(stats_scholtes.iter.X_outer(1,i), stats_scholtes.iter.X_outer(2,i), stats_scholtes.iter.X_outer(1,i+1)-stats_scholtes.iter.X_outer(1,i), stats_scholtes.iter.X_outer(2,i+1)-stats_scholtes.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
    text(stats_scholtes.iter.X_outer(1,i), stats_scholtes.iter.X_outer(2,i),num2str(i));
end
f = @(x,y) x.^2+y.^2-2*(x+y)+2;
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
