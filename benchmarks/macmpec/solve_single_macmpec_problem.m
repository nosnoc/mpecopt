%% load macmpec problems
close all; clear all; clc;
import casadi.*

macmpec_json = dir('macMPEC/*.json');

all_problem_names = {macmpec_json.name};

problame_name = 'qpec-200-2';  
% problame_name =  'ralph1';
% problame_name =  'pack-rig-32';
% problame_name =  'pack-rig-8';
% 
ii_prob = find(contains({macmpec_json.name},problame_name));


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
%% homotopy
settings_homotopy = HomotopySolverOptions();
settings_homotopy.initial_comp_all_zero = 0;
settings_homotopy.homotopy_parameter_steering = "Direct";
% settings_homotopy.comp_tol = 1e-8;
% settings_homotopy.max_iter = 30;

[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
w_opt = full(result_homotopy.x);
f_opt_homotopy = full(result_homotopy.f);
% solver_initalization.w0 = w_opt;
%% mpecopt
fprintf('Problem info, n_w = %d, n_g = %d, n_comp = %d, name = %s\n', length(w),length(g),length(G),name)
solver_settings = MPECOptimizerOptions();
solver_settings.settings_lpec.lpec_solver = 'Highs_casadi';

[results,stats] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_mpecopt = full(results.x);
f_opt_mpecopt = full(results.f);

%% Print result details 
if 1
    fprintf('\n');
    fprintf('Problem info,name = %s, n_w = %d, n_g = %d, n_comp = %d\n',name,length(w),length(g),length(G))
    fprintf('\n-------------------------------------------------------------------------------\n');
    fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf('homotopy \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
    fprintf('mpecopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats.comp_res,stats.n_biactive,stats.cpu_time_total,stats.success,stats.multiplier_based_stationarity)
    fprintf('\n');
    fprintf(' || w_homotopy - w_mpecopt || = %2.2e \n',norm(w_opt-w_opt_mpecopt));
end


