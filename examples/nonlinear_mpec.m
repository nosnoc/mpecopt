clear all
clc
close all
import casadi.*

%%  settings
% violate LICQ
violate_licq = true;
%% Example Nonlinear MPCC
% Variables

% Variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1; x2; x3; x4];

% Nonlinear objective
f = (x1-1)^2 + (x2-1)^2 + (x3-1)^2 + (x4-1)^2 ...
    + 0.5*sin(x1+x2);

% Nonlinear equality (forces pair (x1,x2) biactive possible)
g1 = x1 + x2;       % = 0  → only satisfied if both x1=0 and x2=0
                    % (since x1,x2 ≥ 0 from bounds)

% Nonlinear inequality
g2 = x3^2 + x4^2 - 1;   % ≥ 0  (keeps (x3,x4) on or outside unit circle)


% Complementarity pairs
G = [x1; x3];
H = [x2; x4];

% Bounds
lbx = [0; 0; 0; 0];
ubx = [inf; inf; inf; inf];

% Constraints
g = [g1; g2];
lbg = [0; 0];
ubg = [0; inf];

if violate_licq 
    g = [g; g2];
    lbg = [lbg; 0];
    ubg = [ubg; inf];
end

% Initial guess
x0 = [0.1; 0.1; 1; 0];

mpec = struct('x', x,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);

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
solver_settings.consider_all_complementarities_in_lpec = false;
solver_settings.settings_lpec.lpec_solver = 'Gurobi';
solver_settings.rho_TR_phase_i_init = 10;
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