clear;
clc;
close all;
problem_set_name = 'NLM_small'; % (add large, degenerate, nonlinear, include lifted)
settings.nonlinear_constraints = 1;
settings.copy_ineq = 1; % have two copies of the inequality constrants -> violate licq if active
settings.additional_gradient_term = 1;
settings.include_lifted_vriables = 1;
settings.objective_scale_factor = 1e4;  % 1
settings.s_density = 0.25; % all same or change
settings.objective_functions = {'Linear','Quadratic_psd','Quadratic_ind','Fletcher','Himmelblau','McCormick','Powell','Trigonometric_Nonlinear','Rosenbrock'};
settings.rescale_factor = 10;
settings.n_digits = 6; % significant digits (TODO: can this influence cylcing?)
settings.round_all_data = 0;
settings.bounded_w = 1;
settings.use_normal_distributions_for_feasible_point  = 1;

settings.random_problem_sizes = 1;
settings.NumProbs = 1;  % number of problems of each size
settings.NunProbsRand = 5; % number of problems in range
settings.N_max = 10;
settings.N_min = 2;


mpecs = genNonLinMPECfun(problem_set_name,settings);

%% Define list of solvers to use
solver_names = ["MPEC armin",  'Scholtes direct']; % names of solvers (used for plotting) (strings)
solver_functions = {@mpec_optimizer, @mpec_homotopy_solver};

default_opts1 = MPECOptimizerOptions();
default_opts1.solver_name = solver_names{1};
default_opts1.initialization_strategy = 'RelaxAndProject';
default_opts1.relax_and_project_homotopy_parameter_steering = "Direct";
default_opts1.relax_and_project_iters = 1;
% default_opts1.initialization_strategy = "TakeInitialGuessDirectly";
default_opts1.settings_lpec.lpec_solver = "Gurobi";
default_opts1.debug_mode_on = 0;
default_opts1.verbose_solver = 1;
default_opts1.tol_B_stationarity = 1e-8;
default_opts1.verbose_objective = 0;
default_opts1.stop_if_S_stationary = 0;
% default_opts1.sufficient_decrease_critertion = "Armijo_Linear";
% default_opts1.rho_TR_init = 20; 
% default_opts1.rho_TR_init = 1e6; 
default_opts1.TR_reducing_factor = 0.1;
default_opts1.piece_nlp_strategy  = "BNLP_integer";
default_opts1.rescale_large_objective_gradients = 0;
default_opts1.accept_if_no_objective_improvment  = 0;   

% default_opts1.stop_if_step_rejected= true;
% default_opts1.reset_lpec_objective = false;

scholtes_opts1 = HomotopySolverOptions();
scholtes_opts1.homotopy_parameter_steering = 'Ell_inf';  % 'Ell_inf';
scholtes_opts1.solver_name = solver_names{2};
scholtes_opts1.lift_complementarties = 0;
scholtes_opts1.comp_tol = 1e-9;

opts = {default_opts1, scholtes_opts1}; % list of options to pass to mpecsol (option structs)

%% Create data struct
dstruct = struct;
dstruct.solver_name = [];
dstruct.problem_name = [];
dstruct.success = logical([]);
dstruct.comp_res = [];
dstruct.cpu_time = [];
dstruct.cpu_time_nlp = [];
dstruct.cpu_time_lpec = [];
dstruct.n_biactive = [];
dstruct.n_nlp_total = [];
dstruct.n_lpec_total = [];
dstruct.multiplier_based_stationarity = [];
dstruct.b_stationarity = [];
dstruct.max_iterations_reached = [];
dstruct.problem_infeasible = [];
dstruct.prob_num = [];
dstruct.f = [];
n_problems = length(mpecs);

% find(dtable.max_iterations_reached(dtable.solver_name=='MPEC armin'))
n_mpecs = 1:length(mpecs);
% n_mpecs = [1 15 26]; % max iter
% n_mpecs = 26;
% n_mpecs = [1 11 12 15 18 21]; % infeasible;
% n_mpecs = [10 26 27]; % max iter extended
% n_mpecs = [5 6 14 15 18 21 26]; % no b stat nlp
% n_mpecs = [] 

% n_mpecs = 15;
for ii= [1 2]
    name = solver_names(ii);
    options = opts{ii};
    j = 1;
    % options.relax_and_project_iters
    for jj = n_mpecs
        mpec = mpecs(jj);
        w = mpec.w;
        f = mpec.f_fun(mpec.w);
        g = mpec.g_fun(mpec.w);
        G = mpec.G_fun(mpec.w);
        H = mpec.H_fun(mpec.w);
        w0 = mpec.w0;
        lbw = mpec.lbw;
        ubw = mpec.ubw;
        lbg = mpec.lbg;
        ubg = mpec.ubg;
        mpec_name = mpec.name;
        fprintf('\n\n');
        disp(['solving ' char(mpec.name) ' with solver ' char(name) ' problem: ' num2str(j) '/' num2str(n_problems)]);
        mpec_struct = struct('x',w,'f',f,'g',g,'G',G,'H',H);
        solver_initalization = struct('x0', w0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg, 'ubg',ubg);
        [result,stats] = solver_functions{ii}(mpec_struct,solver_initalization,options);
        dstruct.solver_name = [dstruct.solver_name; string(name)];
        dstruct.problem_name = [dstruct.problem_name; string(mpec_name)];
        dstruct.success = [dstruct.success; logical(stats.success)];
        dstruct.comp_res = [dstruct.comp_res; stats.comp_res];
        dstruct.n_nlp_total = [dstruct.n_nlp_total; stats.n_nlp_total];
        dstruct.n_lpec_total = [dstruct.n_lpec_total; stats.n_lpec_total];
        dstruct.cpu_time = [dstruct.cpu_time; stats.cpu_time_total];
        dstruct.cpu_time_nlp = [dstruct.cpu_time_nlp; stats.cpu_time_nlp];
        dstruct.cpu_time_lpec = [dstruct.cpu_time_lpec; stats.cpu_time_lpec];
        dstruct.n_biactive = [dstruct.n_biactive; stats.n_biactive];
        dstruct.max_iterations_reached = [dstruct.max_iterations_reached ; stats.max_iterations_reached];
        dstruct.problem_infeasible = [dstruct.problem_infeasible; stats.problem_infeasible];
        dstruct.prob_num = [dstruct.prob_num; j];
        dstruct.f = [dstruct.f; result.f];
        dstruct.b_stationarity = [dstruct.b_stationarity; stats.b_stationarity];
        dstruct.multiplier_based_stationarity  = [dstruct.multiplier_based_stationarity; stats.multiplier_based_stationarity];
        j = j+1;
    end
end

%% Do plotting
dtable = struct2table(dstruct);
filename = 'NonlineMPEC_small';

plot_settings.success_fail_statistics = 1;
plot_settings.n_biactive = 0;
plot_settings.lpecs_solved = 0;
plot_settings.nlps_solved = 0;
plot_settings.bar_timing_plots = 1;
plot_settings.objective = 0;
plot_settings.absolute = 1;
plot_settings.relative = 1;
plot_settings.stationary_points = 1;
plot_settings.b_stationarity = 1;
plot_settings.objective_rescaled = 1;
plot_settings.b_stationarty_as_success_criterion = 1;

plot_mpec_benchmark_result(dtable,filename,plot_settings)
stationarity_tab  = create_comparison_table(dtable,{'f','multiplier_based_stationarity'});
stationarity_tab  = create_comparison_table(dtable,{'f','success'});
stationarity_tab  = create_comparison_table(dtable,{'f'});
