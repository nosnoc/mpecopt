% Kintro example  
clear all
clc
close all
import casadi.*

%% Example  MPCC
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x5 = SX.sym('x5');
x6 = SX.sym('x6');
x7 = SX.sym('x7');
x8 = SX.sym('x8');

w = [x1;x2;x3;x4;x5;x6;x7;x8];

p = SX.sym('p');
f = (x1-5)^2+(2*x2+1)^2;
g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
    3*x1-x2-3-x6;...
    -x1+0.5*x2+4-x7;...
    -x1-x2+7-x8];
G = [x6;x7;x8];
H = [x3;x4;x5];
x0 = zeros(8,1);
lbw = zeros(8,1);
ubw = inf*ones(8,1);

lbg = zeros(4,1);
ubg = zeros(4,1);

mpec = struct('x', w, 'f', f, 'g', g,'p',p,'G',G ,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg, 'ubg',ubg,'p0',1);
settings = HomotopySolverOptions();
[sol_homotopy,stats_homotopy] =  mpec_homotopy_solver(mpec,solver_initalization,settings);
x_opt_homotopy = full(sol_homotopy.x);

g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
    3*x1-x2-3-x6;...
    -x1+0.5*x2+4-x7;...
    -x1-x2+7-x8;
    ];
lbg = zeros(4,1);
ubg = zeros(4,1);

%% Homotopy solver
settings_homotopy = HomotopySolverOptions();
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);

%% MINLP solver
settings_minlp = MINLPSolverOptions();
[result_minlp,stats_minlp] = mpec_minlp_solver(mpec,solver_initalization,settings_minlp);
f_opt_minlp = full(result_minlp.f);
w_opt_minlp = full(result_minlp.x);

%% MPECopt solver
solver_settings = mpecopt.Options();
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
solver_settings.settings_lpec.lpec_solver = 'Highs_casadi';
% solver_settings.initialization_strategy = "FeasibilityEll1General";
solver_settings.rho_TR_phase_i_init = 10;
solver_settings.tol_active = 1e-6;
solver = mpecopt.Solver(mpec, solver_settings);
[result_mpecopt,stats_mpecopt] = solver.solve(solver_initalization);
w_opt_mpecopt = full(result_mpecopt.x);
f_opt_mpecopt = full(result_mpecopt.f);

%% Results comparison
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Success\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg     \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('MINLP \t\t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_minlp,stats_minlp.comp_res,stats_minlp.n_biactive,stats_minlp.cpu_time_total,stats_minlp.success,stats_minlp.multiplier_based_stationarity)
fprintf('MPECopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats_mpecopt.comp_res,stats_mpecopt.n_biactive,stats_mpecopt.cpu_time_total,stats_mpecopt.success,stats_mpecopt.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf('||w_reg - w_mpec|| = %2.2e \n',norm(w_opt_homotopy-w_opt_mpecopt));
fprintf('||w_minlp - w_mpec|| = %2.2e \n',norm(w_opt_minlp-w_opt_mpecopt));
fprintf('Solution: (%2.2f,%2.2f) \n',w_opt_mpecopt(1),w_opt_mpecopt(2));

