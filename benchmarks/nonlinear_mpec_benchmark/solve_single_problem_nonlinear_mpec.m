clear; clc; close all;
problem_set_name = 'nonlinear_mpec';
%% Generate test set
settings.nonlinear_eq_constraints = 1;
settings.nonlinear_ineq_constraints = 1;
settings.s_nonlinear_eq = 0.5;
settings.s_nonlinear_ineq = 0.5;
settings.include_lifted_variables_in_objective = 1;
settings.copy_ineq = 1;  %have two copies of the inequality constrants -> violate licq if active
settings.s_ineq_copy = 0.5;
settings.copy_eq = 0;  %have two copies of the inequality constrants -> violate licq if active
settings.s_eq_copy = 0.25;

settings.inequality_constraints = 1;
settings.inequality_constraints_coupling_terms = 1; % if 0, then in Ax+By >=f, B = 0;

settings.objective_functions = {'Quadratic_psd','Quadratic_ind',... 
    'Fletcher','Himmelblau','McCormick',...
    'Powell','Trigonometric','Rosenbrock',...
    'Raydan1','Raydan2',...
    'Diagonal3','Diagonal4','Diagonal5',...
    'Extended_Trigiaonal','Three_Exponential_Terms',...
    'Generalized_PSC1','Extended_PSC1',...
    'Fletcvb3','Bdqrtic','Tridia',...
    'EG2','Edensch','Indef',...
    'Cube','Bdexp','Genhumps',...
    'Arwhead', 'Quartc', 'Cosine', 'Sine'};

settings.rescale_factor = 1;
settings.round_all_data = 0;
settings.n_digits = 6;
settings.bounded_w = 1;
settings.use_normal_distributions_for_feasible_point = 1;
settings.alpha_lin = 10; % f + alpha_lin*[c;d]'*[x;y]'
settings.n_non_zero_E = 2000;
settings.symmetric_psd_M_matrix = 0;

settings.s_density_A_B = 0.1; % all same or change
settings.s_density_M = 0.5;
settings.variable_density = 1;
settings.range_s_density = [0.1 0.3];
% settings.range_s_density = [0.1 0.2];

% settings.s_density_M = 0.1;
% settings.range_s_density = [0.1 0.8];

settings.random_problem_sizes = 1;
% settings.n_ineq_ub = 1.5; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
% settings.n_ineq_lb = 0.5;

settings.n_ineq_ub = 1; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
settings.n_ineq_lb = 0.1;
dimensions.N_rand_prob = 5; % number of problems per objective
dimensions.n_x_max = 100;
dimensions.n_x_min = 10;

settings.n_ineq_ub = 2; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
dimensions.N_rand_prob = 4; % number of problems per objective
dimensions.n_x_max = 50;
dimensions.n_x_min = 10;

% settings.n_ineq_ub = 2; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
% dimensions.N_rand_prob = 1; % number of problems per objective
% dimensions.n_x_max = 30;
% dimensions.n_x_min = 10;

dimensions.n_fraction_of_x = 1; % n_y = round(n_x/n_fraction_of_x)

% example for nonradnom problem sizes
dimensions.n_x_vec = [10,20,50,75,100,100];
dimensions.n_y_vec = [10,10,25,40,50,70];
dimensions.n_ineq_vec = [10,10,20,50,40,250];

mpecs = generate_nonlinear_mpec_problem_set(problem_set_name,settings,dimensions);
length(mpecs)
% pack-rig-32.nl % scholtes

%% mpecopt solvers
   % 2
   

ii_prob = 85; % 10    41    48     85    87
close all;

mpec = mpecs(ii_prob);
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

fprintf('Problem info:name = %s, n_w = %d, n_g = %d, n_comp = %d\n',mpec_name, length(w),length(g),length(G))
mpec_struct = struct('x',w,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', w0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg,'ubg',ubg);
% homotopy
settings_homotopy = HomotopySolverOptions();
settings_homotopy.homotopy_parameter_steering ="Direct";
% settings_homotopy.comp_tol = 1e-12;
% settings_homotopy.max_iter = 25;
% settings_homotopy.lift_complementarities = 1;
% settings_homotopy.lift_complementarities_full = 0;
% settings_homotopy.check_B_stationarity = 1;
% settings_homotopy.comp_res_bilinear = 1;
% settings_homotopy.homotopy_parameter_steering = "Ell_1";
settings_homotopy.plot_mpec_multipliers = true;
% settings_homotopy.comp_tol = 1e-12;
% settings_homotopy.tol_active = 1e-12;
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
w_opt_homotopy = full(result_homotopy.x);
f_opt_homotopy = full(result_homotopy.f);
% solver_initalization.w0 = w_opt_homotopy;
% mpecopt
fprintf('Problem info, n_w = %d, n_g = %d, n_comp = %d, name = %s\n', length(w),length(g),length(G),mpec_name)
solver_settings = MPECOptimizerOptions();
solver_settings.plot_mpec_multipliers = true;
% solver_settings.rho_TR_phase_i_init = 1e1;
solver_settings.stop_if_S_stationary = 0;
% solver_settings.tol_active = 1e-12;
solver_settings.initialization_strategy = "FeasibilityEll1General";
[results,stats] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(results.x);
f_opt_active_set = full(results.f);
% stats.success_phase_i
% stats.success
%
if 1
    fprintf('\n');
    fprintf('Problem info,name = %s, n_w = %d, n_g = %d, n_comp = %d\n',mpec_name,length(w),length(g),length(G))
    fprintf('\n-------------------------------------------------------------------------------\n');
    fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Sucess\t Stat. type\n')
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf('homotopy \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
    fprintf('Active Set \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_active_set,stats.comp_res,stats.n_biactive,stats.cpu_time_total,stats.success,stats.multiplier_based_stationarity)
    fprintf('\n');
    fprintf(' || w_homotopy - w_active_set || = %2.2e \n',norm(w_opt_homotopy-w_opt_active_set));
    fprintf(' || f_homotopy - f_active_set || = %2.2e \n',norm(f_opt_homotopy-f_opt_active_set));
end


%% 
figure
plot(w_opt_active_set-w_opt_homotopy);

