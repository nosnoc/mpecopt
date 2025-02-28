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

solver_settings = mpecopt.Options();
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
solver_settings.initialization_strategy = "FeasibilityEll1General";


 % [sol_active_set,stats_active_set]  = mpec_optimizer(mpec, solver_initalization, solver_settings);
solver = Mpecopt(mpec, solver_settings);
[sol_active_set,stats_active_set] = solver.solve(solver_initalization);


x_opt_active_set = full(sol_active_set.x);
f_opt_homotopy = full(sol_homotopy.f);
f_opt_mpec_opt = sol_active_set.f;


fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Scholtes \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpec_opt,stats_active_set.comp_res,stats_active_set.n_biactive,stats_active_set.cpu_time_total,stats_active_set.success,stats_active_set.multiplier_based_stationarity)
fprintf('\n');
fprintf(' || x_reg- x_active_set || = %2.2e \n',norm(x_opt_homotopy-x_opt_active_set));


