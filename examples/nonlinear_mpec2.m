clear all
clc
close all
import casadi.*

%%  settings
% violate LICQ
violate_licq = true;
import casadi.*

% Variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1; x2; x3; x4];

% Nonlinear objective
f = 5*(x1-1)^2 + (x2-1)^2 + (x3-1)^4 + (x4-1)^4 +10*(x3+x4); 

% Nonlinear equality constraint
g1 = x1+x2+x3+x4;  % = 0

% Nonlinear inequality constraint
g2 = (x1-0.5)^2 + (x3-0.5)^2 - 0.25; % >= 0, keeps vars near origin

% Complementarity pairs
G = [x1; x3];
H = [x2; x4];

% Bounds
lbx = [0; 0; 0; 0];
ubx = [inf; inf; inf; inf];

% Collect constraints
g = [g1; g2];
lbg = [2; 0];     % eq:0, ineq: >=0
ubg = [2; inf];

% Initial guess
x0 = [0.1; 0.1; 0.1; 0.1];

if violate_licq 
    g = [g; g2];
    lbg = [lbg; 0];
    ubg = [ubg; inf];
end

mpec = struct('x', x,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);

%% Reg homotopy
settings_homotopy = HomotopySolverOptions();
settings_homotopy.homotopy_parameter_steering = 'Direct';
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);

%% MINLP reformulation
settings_minlp = MINLPSolverOptions();
[result_minlp,stats_minlp] =  mpec_minlp_solver(mpec,solver_initalization,settings_minlp);
f_opt_minlp = full(result_minlp.f);
w_opt_minlp = full(result_minlp.x);

%%  MpecOpt
solver_settings = mpecopt.Options();
solver_settings.consider_all_complementarities_in_lpec = false;
solver_settings.settings_lpec.lpec_solver = 'Gurobi';
solver_settings.rho_TR_phase_i_init = 10;
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
% solver_settings.initialization_strategy = "FeasibilityEllInfGeneral";

solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% Create solver object
solver = mpecopt.Solver(mpec, solver_settings);
% solve problem
[result_mpecopt,stats_mpecopt] = solver.solve(solver_initalization);

w_opt_mpecopt = full(result_mpecopt.x);
f_opt_mpecopt = full(result_mpecopt.f);

%% print and compare results
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Success\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg     \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('MINLP \t\t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_minlp,stats_minlp.comp_res,stats_minlp.n_biactive,stats_minlp.cpu_time_total,stats_minlp.success,stats_minlp.multiplier_based_stationarity)
fprintf('MPECopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats_mpecopt.comp_res,stats_mpecopt.n_biactive,stats_mpecopt.cpu_time_solvers,stats_mpecopt.success,stats_mpecopt.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf('||w_reg - w_mpec|| = %2.2e \n',norm(w_opt_homotopy-w_opt_mpecopt));
fprintf('||w_minlp - w_mpec|| = %2.2e \n',norm(w_opt_minlp-w_opt_mpecopt));
fprintf('solution is (%2.2f,%2.2f,%2.2f,%2.2f) \n',w_opt_mpecopt(1),w_opt_mpecopt(2),w_opt_mpecopt(3),w_opt_mpecopt(4));
