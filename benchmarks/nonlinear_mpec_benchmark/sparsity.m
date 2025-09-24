import casadi.*

% Define variables
clear; close all;
problem_set_name = 'nonlinear_mpec';
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
density_vec = [];
settings.objective_functions_all = {'Fletcher','McCormick',...
    'Powell','Rosenbrock',...
    'Raydan1','Raydan2',...
    'Diagonal4','Diagonal5',...
    'Extended_Trigiaonal','Three_Exponential_Terms',...
    'Generalized_PSC1',...
    'Fletcvb3','Bdqrtic','Tridia',...
    'EG2','Edensch','Indef',...
    'Cube','Bdexp','Genhumps',...
    'Arwhead', 'Quartc', 'Cosine',...
    'NCVXQP6','DIXCHLNV',...
    };

for jj = 1:27
settings.objective_functions = {settings.objective_functions_all{jj}};

% settings.objective_functions = {'SCURLY30'};



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


settings.nnz_bounded_by_dim = 1;
settings.inv_cond_num = 1e0;
settings.nnz_factor = 1;

settings.adaptive_density_bounds = 0; % to account for very larg problems
settings.variable_density = 1;
settings.range_s_density = [0.01 0.05];
settings.random_problem_sizes = 1;

settings.n_ineq_ub = 2; % n_ineq = n_ineq_ub*n_x % max relative number of ineq w.r.t x
settings.n_ineq_lb = 0.5;

settings.n_fraction_of_x = 0.5;

dimensions.N_rand_prob = 1; % number of problems per objective
dimensions.n_x_max = 300;
dimensions.n_x_min = 300;



dimensions.n_fraction_of_x = 0.5; % n_y = round(n_x/n_fraction_of_x)
mpecs = generate_nonlinear_mpec_problem_set(problem_set_name,settings,dimensions);


w = mpecs.w;
f = mpecs.f;
% Compute Hessian
H = hessian(f, w);
H_func = casadi.Function('H', {w}, {H});
n_obj = length(w);
% Evaluate at test point
w_test = 2.5 * rand(n_obj, 1);
H_val = full(H_func(w_test));
density = nnz(H_val)/length(H_val(:))
density_vec = [density_vec,density];
% Plot sparsity pattern
% figure;
% spy(H_val);
% title('Hessian Sparsity Pattern');
end