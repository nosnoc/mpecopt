clear all
clc
close all
import casadi.*

%%  settings
% Simple nonlinear MPEC test problem
import casadi.*

% Variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];

% Simple quadratic objective
f = (x1-1)^4 + (x2-2)^4;

% Simple equality constraint
g1 = x1 + x2 - 1; % = 0

% Complementarity pairs
G = x1;
H = x2;

% Bounds
lbx = [0; 0];
ubx = [inf; inf];

% Constraints
g = g1;
lbg = 0;
ubg = 0;

% Initial guess
x0 = [0.5; 0.5];

mpec = struct('x', x, 'f', f, 'g', g, 'G', G, 'H', H);
solver_initialization = struct('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);

%% Homotopy solver
settings_homotopy = HomotopySolverOptions();
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initialization,settings_homotopy);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);

%% MINLP solver
settings_minlp = MINLPSolverOptions();
[result_minlp,stats_minlp] = mpec_minlp_solver(mpec,solver_initialization,settings_minlp);
f_opt_minlp = full(result_minlp.f);
w_opt_minlp = full(result_minlp.x);

%% MPECopt solver
solver_settings = mpecopt.Options();
solver_settings.consider_all_complementarities_in_lpec = false;
solver_settings.settings_lpec.lpec_solver = 'Gurobi';
solver = mpecopt.Solver(mpec, solver_settings);
[result_mpecopt,stats_mpecopt] = solver.solve(solver_initialization);
w_opt_mpecopt = full(result_mpecopt.x);
f_opt_mpecopt = full(result_mpecopt.f);


%% Results comparison
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Success\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg     \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('MINLP \t\t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_minlp,stats_minlp.comp_res,stats_minlp.n_biactive,stats_minlp.cpu_time_total,stats_minlp.success,stats_minlp.multiplier_based_stationarity)
fprintf('MPECopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats_mpecopt.comp_res,stats_mpecopt.n_biactive,stats_mpecopt.cpu_time_solvers,stats_mpecopt.success,stats_mpecopt.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf('||w_reg - w_mpec|| = %2.2e \n',norm(w_opt_homotopy-w_opt_mpecopt));
fprintf('||w_minlp - w_mpec|| = %2.2e \n',norm(w_opt_minlp-w_opt_mpecopt));
fprintf('Solution: (%2.2f,%2.2f) \n',w_opt_mpecopt(1),w_opt_mpecopt(2));