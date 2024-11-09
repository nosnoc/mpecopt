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
[sol,stats] =  mpec_homotopy_solver(mpec,solver_initalization,settings);
w_opt = full(sol.x);

g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
    3*x1-x2-3-x6;...
    -x1+0.5*x2+4-x7;...
    -x1-x2+7-x8;
    ];
lbg = zeros(4,1);
ubg = zeros(4,1);

solver_settings = MPECOptimizerOptions();
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
solver_settings.initalization_strategy = "FeasibilityEll1General";


solution = mpec_optimizer(mpec, solver_initalization, solver_settings);
x_opt = full(solution.x);
f_opt_ipopt = full(sol.f)
f_opt_mpec_opt = solution.f


