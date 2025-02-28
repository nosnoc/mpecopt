%% Cardinality optimization

% solve random cardinality optimization problem, using the foromulation
% from: 
% Complementarity Formulations of \ell_0-norm Optimization Problems
% Mingbin Feng, John E. Mitchell, Jong-Shi Pang, Xin Shen, Andreas W achter

clear all; clc; close all
import casadi.*

%% Problem formulation
n = 300; % number of variabples
m = 150;
A = rand(m,n);
b = rand(m,1);
c = 1-2*rand(n,1);
Q = rand(n,n);
e = ones(n,1);
gamma = 1e3;

x = SX.sym('x',n);
x_p = SX.sym('x_p',n);
x_n = SX.sym('x_n',n);
xi = SX.sym('xi',n);
w = [x;x_p;x_n;xi];


f = 0.5*x'*Q*x+c'*x+gamma*e'*(e-xi);

% f = e'*(e-xi);

g = [A*x-b; x-(x_p-x_n)];
x0 = ones(4*n,1);
G = [xi;x_p];
H = [x_p+x_n;x_n];
lbx = -inf(4*n,1);
ubx = 1e2*ones(4*n,1);
ubx(end-n:end) = 1;
lbg = [zeros(m,1);zeros(n,1)];
ubg = [inf(m,1);zeros(n,1)];

mpec = struct('x', w,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
settings = HomotopySolverOptions();
settings.homotopy_parameter_steering = "Direct";
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);
%%  Settings
solver_settings = mpecopt.Options();
% solver_settings.settings_lpec.lpec_solver = "Highs_casadi";
% solver_settings.settings_lpec.lpec_solver = "Gurobi";
solver_settings.rho_TR_phase_i_init = 1e-1;
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
% solver_settings.initalization_strategy = "TakeInitialGuessDirectly";

solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
%%
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('homotopy \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf(' || w_homotopy - w_active_set || = %2.2e \n',norm(w_opt_homotopy-w_opt_active_set));

x_opt_active_set = w_opt_active_set(1:n);
x_opt_homotopy = w_opt_homotopy(1:n);
%%
cardinality_active_set = sum(heaviside(abs(x_opt_active_set)-1e-3))
cardinality_homotopy  = sum(heaviside(abs(x_opt_homotopy)-1e-3))
