clear all; clc; close all;
import casadi.*
%% Kirches 2022, example, origin is M stationary, (0,1) is the global optimum

x1 = SX.sym('x1');
x2 = SX.sym('x2');

x = [x1;x2];

% parameters
p = SX.sym('p');
p0 = 0.5;

f = x1^3-x2+p*x2^2;
G = x1;
H = x2;

x0 = [0;5];
x0 = [3;5];
x0 = [3;0];
lbx = [0;0];
ubx = [inf;inf];

g = [x1^2+x2^2];
lbg = [-inf];
ubg = [15];

g = [ ];
lbg = [ ];
ubg = [ ];

mpec = struct('x', x, 'f', f, 'g', g,'G',G,'H',H,'p',p);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p0',p0);
% Scholtes
settings = HomotopySolverOptions();
[result_scholtes,stats_scholtes] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_scholtes = full(result_scholtes.f);
x_opt_scholtes = full(result_scholtes.x);

% x0 = w_opt_scholtes;
solver_settings = mpecopt.Options();
% solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
solver_settings.consider_all_complementarities_in_lpec = true;
solver_settings.tol_B_stationarity = 1e-8;
solver_settings.rho_TR_phase_ii_init = 1e-2;
% solver_initalization.x0 = x0;

solver = mpecopt.Solver(mpec, solver_settings);
[result_active_set,stats_active_set] = solver.solve(solver_initalization);


x_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);


fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Scholtes \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_scholtes,stats_scholtes.comp_res,stats_scholtes.n_biactive,stats_scholtes.cpu_time_total,stats_scholtes.success,stats_scholtes.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('\n');
fprintf(' || x_scholtes - x_active_set || = %2.2e \n',norm(x_opt_scholtes-x_opt_active_set));

