clear; clc; close all;
problem_set_name = 'nonlinear_mpec2'; % (add large, degenerate, nonlinear, include lifted)

%% file names
filename = 'random_lpecs2'; % name for figures and excel table
results_name = ['results/' filename '_' datestr(datetime("today"))]; % name for matlab .dat with results
%% Generate test set
settings.casadi_variable_type = 'SX';
settings.nonlinear_eq_constraints = 1;
settings.nonlinear_ineq_constraints = 1;
settings.s_nonlinear_eq = 0.01;
settings.s_nonlinear_ineq = 1/4;
settings.include_lifted_variables_in_objective = 0;
settings.copy_ineq = 1;  %have two copies of the inequality constrants -> violate licq if active
settings.s_ineq_copy = 1/4;
settings.copy_eq = 0;  
settings.s_eq_copy = 0.0;

settings.inequality_constraints = 1;
settings.inequality_constraints_coupling_terms = 1; % if 0, then in Ax+By >=f, B = 0;

% Trigonometric
settings.objective_functions = {'Quadratic_psd','Quadratic_ind',... 
    'Fletcher','Himmelblau','McCormick',...
    'Powell','Rosenbrock',...
    'Raydan1','Raydan2',...
    'Diagonal4','Diagonal5',...
    'Extended_Trigiaonal','Three_Exponential_Terms',...
    'Generalized_PSC1','Extended_PSC1',...
    'Fletcvb3','Bdqrtic','Tridia',...
    'EG2','Edensch','Indef',...
    'Cube','Bdexp','Genhumps',...
    'Arwhead', 'Quartc', 'Cosine', 'Sine',...
    'CURLY20', 'SCURLY30', 'FLETBV3M',...
    'NCVXQP6','DIXCHLNV','EDENSCH',...
    };


settings.rescale_factor = 1;
settings.round_all_data = 1;
settings.n_digits = 4;
settings.bounded_w = 1;
settings.use_normal_distributions_for_feasible_point = 1;
settings.alpha_lin = 10; % f + alpha_lin*[c;d]'*[x;y]'
settings.n_non_zero_E = 2000;
settings.symmetric_psd_M_matrix = 0;

settings.s_density_A_B = 0.1; % all same or change
settings.s_density_M = 0.5;

settings.s_density_A_B = 0.05; % all same or change
settings.s_density_M = 0.1;

settings.variable_density = 1;
settings.range_s_density = [0.01 0.1];
settings.random_problem_sizes = 1;

settings.n_ineq_ub = 3; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
settings.n_ineq_lb = 0.5;

% Problem size 

dimensions.N_rand_prob = 3; % number of problems per objective
dimensions.n_x_max = 425;
dimensions.n_x_min = 100;


% dimensions.N_rand_prob = 3; % number of problems per objective
% dimensions.n_x_max = 50;
% dimensions.n_x_min = 10;

% Large
% dimensions.n_x_max = 1600;
% dimensions.n_x_min = 100;


dimensions.n_fraction_of_x = 0.5; % n_y = round(n_x/n_fraction_of_x)

% example for nonradnom problem sizes
% dimensions.n_x_vec = [10,20,50,75,100,100];
% dimensions.n_y_vec = [10,10,25,40,50,70];
% dimensions.n_ineq_vec = [10,10,20,50,40,250];

mpecs = generate_nonlinear_mpec_problem_set(problem_set_name,settings,dimensions);
length(mpecs)


%% Solver settings
%% Define list of solvers to use
solver_names  = ["Gurobi", "Gurobi-early-I", "Gurobi-early-I-II", "Highs", "Highs-early"];

solver_functions = {@mpec_optimizer,@mpec_optimizer,@mpec_optimizer,@mpec_optimizer,@mpec_optimizer};

opts1 = mpecopt.Options();
opts1.solver_name = solver_names{1};
opts1.settings_lpec.lpec_solver = "Gurobi";
opts1.relax_and_project_homotopy_parameter_steering = "Direct";
opts1.settings_lpec.stop_lpec_at_feasible = false;
opts1.use_one_nlp_solver = true;


opts2 = mpecopt.Options();
opts2.solver_name = solver_names{2};
opts2.settings_lpec.lpec_solver = "Gurobi";
opts2.relax_and_project_homotopy_parameter_steering = "Direct";
opts2.settings_lpec.stop_lpec_at_feasible = true;
opts2.use_one_nlp_solver = true;

opts3 = mpecopt.Options();
opts3.solver_name = solver_names{3};
opts3.settings_lpec.lpec_solver = "Gurobi";
opts3.relax_and_project_homotopy_parameter_steering = "Direct";
opts3.settings_lpec.stop_lpec_at_feasible = false;
opts3.settings_lpec.stop_lpec_at_descent = false;
opts3.use_one_nlp_solver = true;


opts4 = mpecopt.Options();
opts4.solver_name = solver_names{4};
opts4.settings_lpec.lpec_solver = "Highs";
opts4.relax_and_project_homotopy_parameter_steering = "Direct";
opts4.settings_lpec.stop_lpec_at_feasible = false;
opts4.use_one_nlp_solver = true;

opts5 = mpecopt.Options();
opts5.solver_name = solver_names{5};
opts5.settings_lpec.lpec_solver = "Highs";
opts5.relax_and_project_homotopy_parameter_steering = "Direct";
opts5.settings_lpec.stop_lpec_at_feasible = true;
opts5.use_one_nlp_solver = true;

opts = {opts1, opts2, opts3, opts4, opts5}; % list of options to pass to mpecsol (option structs)

%% Create data struct
N_experiments = [1:3];

nonlinear_mpec_benchmark_dtable_loop; % this script runs the experimetns, creates a dtable
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
plot_settings.success_fail_statistics = 0;
plot_settings.n_biactive = 0;
plot_settings.lpecs_solved = 0;
plot_settings.active_set_changes = 0;
plot_settings.lpecs_cpu_time = 0;
plot_settings.bar_timing_plots = 0;
plot_settings.nlps_solved = 0;
plot_settings.max_nlp_cpu_time = 0;
plot_settings.max_nlp_cpu_time_phase_i = 0;
plot_settings.max_nlp_cpu_time_phase_ii = 0;
plot_settings.nlp_cpu_time = 0; % aggegated nlp times phase I and II
plot_settings.nlp_lpec_cpu_comparisson = 0;
plot_settings.objective = 0;
plot_settings.objective_rescaled = 1;
plot_settings.max_lpec_cpu_time = 0;

plot_settings.relative = 1;
plot_settings.relative_phase_i = 0;
plot_settings.relative_phase_ii = 0;
plot_settings.relative_lpec = 0;

plot_settings.absolute = 1;
plot_settings.absolute_phase_i = 0;
plot_settings.absolute_phase_ii = 0;
plot_settings.absolute_lpec = 0;

plot_settings.solved_in_phase_i = 0;

% split cpu time among phases and show in bar plot (just report for mpecopt
% guroi)
plot_settings.nlp_phases_cpu_time = 0;
plot_settings.lpec_phases_cpu_time = 0;

plot_settings.max_lpec_cpu_time = 0;
plot_settings.stationary_points = 0;
plot_settings.b_stationarity = 0;
plot_settings.b_stationarty_as_success_criterion = 0;
plot_settings.plot_only_sucessful = 1;
plot_settings.bar_comparisson_plots = 0;

plot_settings.save_plot = 1;

close all;
plot_mpec_benchmark_result(dtable_plot,filename,plot_settings)
%%
dtable_ii = dtable_plot(dtable_plot.solver_name == solver_names{1},:);
dtable_jj = dtable_plot(dtable_plot.solver_name == solver_names{4},:);
hh = find(dtable_jj.b_stationarity == 0 & dtable_jj.success == 1);
% dtable_ii.f_lpec((hh))
% dtable_jj.f((hh))
% dtable_jj.f_lpec((hh))
% abs((dtable_ii.f((hh))-dtable_jj.f((hh))))./abs(dtable_ii.f((hh)))>1e-3

%%
% figure
% semilogy(abs(dtable_ii.f_lpec)+1e-16,'LineWidth',1.5);
% hold on
% semilogy(abs(dtable_jj.f_lpec)+1e-16,'LineWidth',1.5);
% ylabel('$|\nabla f(x^*)^\top d|$')
% xlabel('Problem instance')
% grid on
% yline(1e-8,'LineWidth',2)