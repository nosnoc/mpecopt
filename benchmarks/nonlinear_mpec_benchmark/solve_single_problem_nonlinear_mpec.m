clear; clc; close all;
problem_set_name = 'nonlinear_mpec';
settings.nonlinear_eq_constraints = 1;
settings.nonlinear_ineq_constraints = 1;
settings.s_nonlinear_eq = 0.05;
settings.s_nonlinear_ineq = 0.05;
settings.include_lifted_variables_in_objective = 0;
settings.copy_ineq = 1;  %have two copies of the inequality constrants -> violate licq if active
settings.s_ineq_copy = 0.2;
settings.copy_eq = 0;  %have two copies of the inequality constrants -> violate licq if active
settings.s_eq_copy = 0.0;

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

settings.objective_functions = {'Himmelblau'};

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
settings.range_s_density = [0.1 0.3];
settings.random_problem_sizes = 1;

settings.n_ineq_ub = 2; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
settings.n_ineq_lb = 0.1;

% Problem size 

dimensions.N_rand_prob = 1; % number of problems per objective
dimensions.n_x_max = 400;
dimensions.n_x_min = 400;

dimensions.n_fraction_of_x = 1; % n_y = round(n_x/n_fraction_of_x)

% example for nonradnom problem sizes
dimensions.n_x_vec = [10,20,50,75,100,100];
dimensions.n_y_vec = [10,10,25,40,50,70];
dimensions.n_ineq_vec = [10,10,20,50,40,250];

mpecs = generate_nonlinear_mpec_problem_set(problem_set_name,settings,dimensions);
length(mpecs)

%% mpecopt solvers
ii_prob = 1; % 10    41    48     85    87
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

%% Set up problem
fprintf('Problem info:name = %s, n_w = %d, n_g = %d, n_comp = %d\n',mpec_name, length(w),length(g),length(G))
mpec_struct = struct('x',w,'f',f,'g',g,'G',G,'H',H);
solver_initalization = struct('x0', w0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg,'ubg',ubg);
%% Homotopy solver
settings_homotopy = HomotopySolverOptions();
settings_homotopy.homotopy_parameter_steering = "Direct";
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings_homotopy);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);

%% MINLP solver
% settings_minlp = MINLPSolverOptions();
% [result_minlp,stats_minlp] = mpec_minlp_solver(mpec,solver_initalization,settings_minlp);
% f_opt_minlp = full(result_minlp.f);
% w_opt_minlp = full(result_minlp.x);

%% MPECopt solver
solver_settings = mpecopt.Options();
solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
solver_settings.settings_lpec.lpec_solver = 'Gurobi';
% solver_settings.rho_TR_phase_i_init = 10;

solver = mpecopt.Solver(mpec, solver_settings);
[result_mpecopt,stats_mpecopt] = solver.solve(solver_initalization);
w_opt_mpecopt = full(result_mpecopt.x);
f_opt_mpecopt = full(result_mpecopt.f);

%% Results comparison
fprintf('\n-------------------------------------------------------------------------------\n');
fprintf('Method \t\t Objective \t comp_res \t n_biactive \t CPU time (s)\t Success\t Stat. type\n')
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Reg     \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_homotopy,stats_homotopy.comp_res,stats_homotopy.n_biactive,stats_homotopy.cpu_time_total,stats_homotopy.success,stats_homotopy.multiplier_based_stationarity)
% fprintf('MINLP \t\t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_minlp,stats_minlp.comp_res,stats_minlp.n_biactive,stats_minlp.cpu_time_total,stats_minlp.success,stats_minlp.multiplier_based_stationarity)
fprintf('MPECopt \t %2.2e \t %2.2e \t\t %d \t\t\t %2.2f \t\t\t\t %d\t %s\n',f_opt_mpecopt,stats_mpecopt.comp_res,stats_mpecopt.n_biactive,stats_mpecopt.cpu_time_total,stats_mpecopt.success,stats_mpecopt.multiplier_based_stationarity)
fprintf('-------------------------------------------------------------------------------\n');
fprintf('||w_reg - w_mpec|| = %2.2e \n',norm(w_opt_homotopy-w_opt_mpecopt));
fprintf('||f_reg - f_mpec|| = %2.2e \n',norm(f_opt_homotopy-f_opt_mpecopt));
% fprintf('||w_minlp - w_mpec|| = %2.2e \n',norm(w_opt_minlp-w_opt_mpecopt));
fprintf('Solution: (%2.2f,%2.2f) \n',w_opt_mpecopt(1),w_opt_mpecopt(2));


