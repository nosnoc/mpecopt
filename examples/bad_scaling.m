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

f = 0.1*((x1-5)^2+(2*x2+1)^2+x4^4+x5^3);
g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
    x1^2+x2^2+x3^2-10
    -x1+0.5*x2+4-x7^2;...
    -x1-x2^2+1-x8^3;...
    x1^2+x2^2+x3^2+x4^2-16;...
    ];
G = [x6;x7;x8];
H = [x3;x4;x5];

x0 = 5*ones(8,1);

lbw = [-2;-2;zeros(6,1)];
ubw = [1;2;3;inf*ones(5,1)];

lbg = [zeros(4,1);-inf;];
ubg = zeros(5,1);

mpec = struct('x', w, 'f', f, 'g', g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
% Scholtes
settings = HomotopySolverOptions();
[result_scholtes,stats_scholtes] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_scholtes = full(result_scholtes.f);
w_opt_scholtes = full(result_scholtes.x);

x0 = w_opt_scholtes;
x0 = [1;1;1;1;0;0;0;0];


solver_settings = MPECOptimizerOptions();
% solver_settings.rho_TR_phase_i_init = 1e-3;
solver_settings.tol_B_stationarity_early_term = 1e-6;
% solver_settings.max_inner_iter = 7;
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);


fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Scholtes \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_scholtes,stats_scholtes.comp_res,stats_scholtes.n_biactive,stats_scholtes.cpu_time_total,stats_scholtes.success,stats_scholtes.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('\n');
fprintf(' || w_scholtes - w_active_set || = %2.2e \n',norm(w_opt_scholtes-w_opt_active_set));

