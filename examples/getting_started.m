clear all
clc
close all
import casadi.*

%% Example  MPCC
% optimization variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
w = [x1;x2;x3];
% parameter
p = SX.sym('p'); 
% objective
f = x1+x2-p*x3;

% inital guess
x0 = ones(3,1);
% bounds
lbw = [0;0;-inf];
ubw = [inf;inf;inf];

% genearl equality and inequality constraints
g = [4*x1-x3;...
     4*x2-x3];
lbg = [0;0];
ubg = [inf;inf];

% complementarity variables
G = x1;
H = x2;

mpec = struct('x', w, 'f', f, 'g', g,'p',p,'G',G ,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg, 'ubg',ubg,'p0',1);

% homotopy' global regularization (reg)
settings = HomotopySolverOptions();
settings.settings_casadi_nlp.ipopt.linear_solver = 'mumps'; % default; 
settings.homotopy_parameter_steering = 'Direct'; % Direct = Scholtes' global regularization, other options: ell_1 and ell_inf penality;
settings.comp_tol = 1e-14; % needs very high accuracy for success
settings.max_iter = 25;
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
w_opt_homotopy = full(result_homotopy.x);
f_opt_homotopy = full(result_homotopy.f);
%
% Pivoting
solver_settings = MPECOptimizerOptions();
% change some settings
solver_settings.settings_lpec.lpec_solver = 'Highs'  ; % 'Gurobi'; for best perfomance;
solver_settings.settings_casadi_nlp.ipopt.linear_solver = 'mumps'; % 'ma27' for better perfomance
solver_settings.rho_TR_phase_i_init = 1e1;
solver_settings.rho_TR_phase_ii_init = 1e-4;

[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);

fprintf('\n');
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg        \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('\n');
fprintf(' || w_homotopy - w_active_set || = %2.2e \n',norm(w_opt_homotopy-w_opt_active_set));
