%% load macmpec problems
clear ; clc; close all;
import casadi.*

%% Color order

results_name = ['results/macmpec_general_' datestr(datetime("today"))]; % name for matlab .dat with results
filename = 'macmpec_general2'; % name for figures and excel table

%% Load macmpec
macmpec_json = dir('macMPEC/*.json');
mpecs = [];
n_problems = length(macmpec_json);

N_not_B_scholtes = [164 166 183 184 133]; % first four low b stat coef
N_biactive = [168 80 83 71 85 73 120 117 105 108];
N_non_S = [84 72 141 123 132 144 147 125 143 134 35 170 181  169];
N_qpecs = [160 161 162 163 164 165 166 167];
N_failed_by_scohltes = [94 95 110 113];
N_not_presolve  = [45 46 47 191];
N_infeasible = [60 66];
N_almost_feasible = [114 68];
N_easy = [34 36 65 87 49 67 69 33 51 53 50 55 57 174 25 40 101 16 102 1 2 3 4]; % bunch of easy problems for more diversity
N_interesting = [N_failed_by_scohltes, N_biactive, N_non_S, N_qpecs, N_not_presolve, N_easy, N_infeasible, N_almost_feasible];
% N_interesting = [N_easy];
% N_interesting = [N_easy, N_failed_by_scohltes, N_biactive];

% for ii = N_interesting
for ii=1:length(macmpec_json)
    fname = fullfile(macmpec_json(ii).folder, macmpec_json(ii).name);
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
    mpecs = [mpecs,mpec];
end

%% Define list of solvers to use

solver_names  = ["MPECopt-Gurobi", "MPECopt-HiGHS",...
                      "MPECopt-Simple", "Reg" , "Pen-$\ell_{\infty}$", "Pen-$\ell_{1}$"];


solver_names  = ["MPECopt-Reg-Gurobi", "MPECopt-Reg-HiGHS",...
                      "MPECopt-$\ell_1$-Gurobi", "Reg" , "Pen-$\ell_{\infty}$", "Pen-$\ell_{1}$"];


solver_functions = {@mpec_optimizer,@mpec_optimizer,@mpec_optimizer,...
        @mpec_homotopy_solver,@mpec_homotopy_solver,@mpec_homotopy_solver};

default_opts1 = MPECOptimizerOptions();
default_opts1.solver_name = solver_names{1};

default_opts2 = MPECOptimizerOptions();
default_opts2.solver_name = solver_names{2};
default_opts2.settings_lpec.lpec_solver = "Highs";
default_opts2.rho_TR_phase_i_init = 1e-3;

default_opts3 = MPECOptimizerOptions();
default_opts3.solver_name = solver_names{3};
default_opts3.relax_and_project_homotopy_parameter_steering = "Ell_1";
% default_opts3.bnlp_projection_strategy = "Simple";

scholtes_opts1 = HomotopySolverOptions();
scholtes_opts1.homotopy_parameter_steering = 'Direct';

scholtes_opts2 = HomotopySolverOptions();
scholtes_opts2.homotopy_parameter_steering = 'Ell_inf';

scholtes_opts3 = HomotopySolverOptions();
scholtes_opts3.homotopy_parameter_steering = 'Ell_1';

opts = {default_opts1, default_opts2, default_opts3, scholtes_opts1, scholtes_opts2, scholtes_opts3}; % list of options to pass to mpecsol (option structs)

%% Create data struct
N_experiments = [1,2,3,4,6];
mpec_benchmark_dtable_loop; % this script runs the experimetns, creates a dtable
%%  Pick which results to plot
dtable = dtable1;
% Modift possibly which solvers are plotted
[~,temp_b] = unique(dtable.solver_name);
solver_name = dtable.solver_name(sort(temp_b));
% solver_tab = table(solver_name);
N_plot = N_experiments;
% N_plot = [1,2];

dtable_plot = [];
for ii = N_plot
    dtable_plot  = [dtable_plot; dtable(dtable.solver_name == solver_names(ii),:)];
end

%% Export to excel 
write_in_separate_sheets = true;
filename_excel = [filename '_details.xlsx'];
if write_in_separate_sheets
    solver_tab = pivot(dtable_plot, Rows= "solver_name");
    for ii = 1:length(solver_tab.solver_name)
        solver_name =solver_tab.solver_name{ii};
        dtable_ii = dtable_plot(dtable_plot.solver_name == solver_name,:);
        writetable(dtable_ii,filename_excel,'Sheet',ii)
    end
else
    writetable(dtable_plot,filename_excel,'Sheet',1)
end

%% Plot results  
plot_settings.success_fail_statistics = 1;
plot_settings.n_biactive = 0;
plot_settings.lpecs_solved = 0;
plot_settings.active_set_changes = 0;
plot_settings.lpecs_cpu_time = 0;
plot_settings.bar_timing_plots = 0;
plot_settings.nlps_solved = 1;
plot_settings.max_nlp_cpu_time = 1;
plot_settings.max_nlp_cpu_time_phase_i = 1;
plot_settings.max_nlp_cpu_time_phase_ii = 1;
plot_settings.nlp_cpu_time = 1; % aggegated nlp times phase I and II
plot_settings.nlp_lpec_cpu_comparisson = 0;
plot_settings.objective = 0;
plot_settings.objective_rescaled = 1;
plot_settings.max_lpec_cpu_time = 0;

plot_settings.relative = 1;
plot_settings.relative_phase_i = 1;
plot_settings.relative_phase_ii = 0;
plot_settings.relative_lpec = 0;

plot_settings.absolute = 1;
plot_settings.absolute_phase_i = 0;
plot_settings.absolute_phase_ii = 0;
plot_settings.absolute_lpec = 0;

plot_settings.solved_in_phase_i = 1;

% split cpu time among phases and show in bar plot (just report for mpecopt
% guroi)
plot_settings.nlp_phases_cpu_time = 0;
plot_settings.lpec_phases_cpu_time = 0;

plot_settings.max_lpec_cpu_time = 0;
plot_settings.stationary_points = 1;
plot_settings.b_stationarity = 0;
plot_settings.b_stationarty_as_success_criterion = 0;
plot_settings.plot_only_sucessful = 1;
plot_settings.bar_comparisson_plots = 0;

plot_settings.save_plot = 1;

close all;
plot_mpec_benchmark_result(dtable_plot,filename,plot_settings)


