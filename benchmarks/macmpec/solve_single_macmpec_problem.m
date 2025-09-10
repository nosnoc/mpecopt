%% load macmpec problems
close all; clear all; clc;
import casadi.*

macmpec_json = dir('macMPEC/*.json');

all_problem_names = {macmpec_json.name};

problame_name = 'qpec-200-3';  
% problame_name =  'ralph1';
% problame_name =  'pack-rig-32';
% problame_name =  'pack-rig-8';
% problame_name = 'gnash15m.nl';
% problame_name = 'pack-rig2p-16.nl';
% problame_name =  'nash1a'; 
% problame_name =  'nash1c'; 
% problame_name =  'tap-09'; 
% problame_name = ' design-cent-31';

% problame_name = 'pack-rig2p-16';

% problame_name = 'pack-rig2p-8';

 % problame_name = 'pack-rig1p-32'; % lot of lpec itters


% Infesaile : gnash15m gnash16m gnash17m gnash18m gnash19m 
% Cyicling: tap-09?
% A-stats: pack-rig2p-16.nl 
% siouxfls1 and siouxfls % bad tolerances?

N_biactive = [168 80 83 71 85 73 120 117 105 108];

ii_prob = find(contains({macmpec_json.name},problame_name));
% ii_prob = N_biactive(10); % A stationarity in reg! looks problematic!
% ii_prob = N_biactive(10);
% ii_prob = 17;


fname = fullfile(macmpec_json(ii_prob).folder, macmpec_json(ii_prob).name);
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
mpec = jsondecode(str);

mpec.w = SX.deserialize(mpec.w);
mpec.f_fun = Function.deserialize(mpec.f_fun);
mpec.g_fun = Function.deserialize(mpec.g_fun);
mpec.G_fun = Function.deserialize(mpec.G_fun);
mpec.H_fun = Function.deserialize(mpec.H_fun);
name = mpec.name

%% create casadi problem
w = mpec.w;
f = mpec.f_fun(mpec.w);
g = mpec.g_fun(mpec.w);
G = mpec.G_fun(mpec.w);
H = mpec.H_fun(mpec.w);

for ii=1:length(H)
    if mpec.lbH(ii) ~= -inf && mpec.ubH(ii) ~= inf
        Hi = H(ii);
        Gi = G(ii);
        H(ii) = Hi - mpec.lbH(ii);
        H = vertcat(H, mpec.ubH(ii) - Hi);
        G = vertcat(G,  -Gi);
    elseif mpec.lbH(ii) == -inf && mpec.ubH(ii) == inf
        error("Something is very wrong")
    elseif mpec.lbH(ii) == -inf
        Hi = H(ii);
        Gi = G(ii);
        H(ii) = mpec.ubH(ii) - Hi;
        G(ii) = -Gi;
    elseif mpec.ubH(ii) == inf
        Hi = H(ii);
        Gi = G(ii);
        H(ii) = Hi - mpec.lbH(ii);
    end
end

% mpecopt solvers
w0 = mpec.w0;
lbw = mpec.lbw;
ubw = mpec.ubw;
lbg = mpec.lbg;
ubg = mpec.ubg;
name = mpec.name;
mpec = struct('x', w, 'f', f, 'g', g,'G',G,'H',H);
solver_initalization = struct('x0', w0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg, 'ubg',ubg);
clc
fprintf('Problem info, n_w = %d, n_g = %d, n_comp = %d, name = %s\n', length(w),length(g),length(G),name)

%% Homotopy solver
settings_homotopy = HomotopySolverOptions();
settings_homotopy.homotopy_parameter_steering = 'Direct';
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);

%% MINLP solver
settings_minlp = MINLPSolverOptions();
settings_minlp.settings_casadi_nlp.bonmin.time_limit = 25;
settings_minlp.settings_casadi_nlp.bonmin.node_limit = 5;
settings_minlp.settings_casadi_nlp.bonmin.solution_limit = 5;
settings_minlp.settings_casadi_nlp.bonmin.max_consecutive_failures = 5;
% [result_minlp,stats_minlp] = mpec_minlp_solver(mpec,solver_initalization,settings_minlp);
% f_opt_minlp = full(result_minlp.f);
% w_opt_minlp = full(result_minlp.x);

%% MPECopt solver
solver_settings = mpecopt.Options();
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
solver_settings.settings_lpec.lpec_solver = "Gurobi";
% solver_settings.settings_casadi_nlp.ipopt.max_iter = 4000;
solver_settings.settings_lpec.stop_lpec_at_feasible = true;
solver_settings.settings_lpec.stop_lpec_at_descent = true;
% solver_settings.settings_casadi_nlp.ipopt.fixed_variable_treatment = 'relax_bounds';
% solver_settings.initialization_strategy = "FeasibilityEll1General";
% solver_settings.rho_TR_phase_ii_init = 1e-4;
% solver_settings.consider_all_complementarities_in_lpec = false;
% solver_settings.tol_active = 1e-6;
solver_settings.use_one_nlp_solver = true;
solver_settings.compute_tnlp_stationary_point = false;
tic
solver = mpecopt.Solver(mpec, solver_settings);
toc

[result_mpecopt,stats_mpecopt] = solver.solve(solver_initalization);
stats_mpecopt.cpu_time_total

% [solution,stats] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_mpecopt = full(result_mpecopt.x);
f_opt_mpecopt = full(result_mpecopt.f);

stats_mpecopt.iter.nodecount_phase_i
stats_mpecopt.iter.nodecount_phase_ii
stats_mpecopt.iter.cpu_time_lpec_phase_i_iter
stats_mpecopt.iter.cpu_time_lpec_phase_ii_iter

%% Results comparison
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Success\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg     \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
% fprintf('MINLP \t\t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_minlp,stats_minlp.comp_res,stats_minlp.n_biactive,stats_minlp.cpu_time_total,stats_minlp.success,stats_minlp.multiplier_based_stationarity)
fprintf('MPECopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats_mpecopt.comp_res,stats_mpecopt.n_biactive,stats_mpecopt.cpu_time_total,stats_mpecopt.success,stats_mpecopt.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf('||w_reg - w_mpec|| = %2.2e \n',norm(w_opt_homotopy-w_opt_mpecopt));
% fprintf('||w_minlp - w_mpec|| = %2.2e \n',norm(w_opt_minlp-w_opt_mpecopt));


% 
% stats_mpecopt.iter.itercount_phase_i
% stats_mpecopt.iter.itercount_phase_ii