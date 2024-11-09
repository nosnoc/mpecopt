clear all; clc; 
close all;
import casadi.*
%% Kirches 2022 Example from Section 5.3
solve_as_nlp = 1;
use_active_set_method = 1;
x = SX.sym('x',8);
x_0 = x(1:4);
x_1 = x(5:6);
x_2 = x(7:8);
rho = 2;
lambda1= 3.9375;
lambda2 = -6.5;
lambda3 = -0.25;
lambda4 = 2.5;

if solve_as_nlp
    g = [-34+2*x_0(3)+8/3*x_0(4)+x_2(1);
        (-24.25+1.25*x_0(3)+2*x_0(4)+x_2(2));
        (x_1(1)+x_0(2)+x_0(3)-15);
        x_1(2)+x_0(1)-x_0(4)-15];
    lbg = zeros(4,1);
    ubg = inf*ones(4,1);
    f = 0.5*((x_0(1)-x_0(3))^2+(x_0(2)-x_0(4))^2);
else
    f = 0.5*((x_0(1)-x_0(3))^2+(x_0(2)-x_0(4))^2)...
        +lambda1*(-34+2*x_0(3)+8/3*x_0(4)+x_2(1))...
        -lambda2*(-24.25+1.25*x_0(3)+2*x_0(4)+x_2(2))...
        -lambda3*(x_1(1)+x_0(2)+x_0(3)-15)...
        +lambda4*(x_1(2)+x_0(1)-x_0(4)-15)...
        +0.5*rho*( (-34+2*x_0(3)+8/3*x_0(4)+x_2(1))^2 + (-24.25+1.25*x_0(3)+2*x_0(4)+x_2(2))^2+ (x_1(1)+x_0(2)+x_0(3)-15)^2+ (x_1(2)+x_0(1)-x_0(4)-15)^2 );
    g = [];
    lbg = [];
    ubg = [];
end

G = x_1;
H = x_2;

x0 = zeros(8,1);
% x0 = 50*rand(8,1);

% x0 = [  11.6667
%     3.2083
%    11.6667
%     3.2083
%     0.0000
%     5.2917
%     0.1424
%    -0.0000];

lbx = [-inf*ones(4,1);0*ones(4,1)];
ubx = inf*ones(8,1);

mpec = struct('x', x, 'f', f, 'g', g,'G',G ,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx,'ubx',ubx,'lbg',lbg, 'ubg',ubg,'p0',1);

% Scholtes
settings = HomotopySolverOptions();
[result_scholtes,stats_scholtes] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_scholtes = full(result_scholtes.f);
w_opt_scholtes = full(result_scholtes.x);

% x0 = w_opt_scholtes;
solver_settings = MPECOptimizerOptions();
solver_settings.initalization_strategy = "TakeInitialGuessDirectly";
solver_settings.consider_all_complementarities_in_lpec = true;
solver_settings.settings_lpec.lpec_solver = "Projected_Gradient";
% solver_initalization.x0 = x0;

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

