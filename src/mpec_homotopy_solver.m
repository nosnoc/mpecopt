function [solution,stats] = mpec_homotopy_solver(mpec,solver_initalization,settings)
import casadi.*
try
    x = mpec.x;
catch
    x = mpec.w;
end
f = mpec.f;
g = mpec.g;
G = mpec.G;
H = mpec.H;
x_k = solver_initalization.x0;
lbx = solver_initalization.lbx;
ubx = solver_initalization.ubx;
lbg = solver_initalization.lbg;
ubg = solver_initalization.ubg;

if isfield(mpec,'p')
    p = mpec.p;
    p0 = solver_initalization.p0;
else
    p = [];
    p0 = [];
end

%% homotopy parameters
sigma = SX.sym('sigma',1);
ii = 1;
comp_res_ii = 1e3;
p = [p;sigma];
p0 = [p0;settings.sigma0];
%% comp functions
n_primal_non_lifted = length(x);
n_comp = size(G,1);
G_fun = Function('G_fun',{x,p},{G});
H_fun = Function('G_fun',{x,p},{H});

G_copy = G_fun(x,p);
H_copy = H_fun(x,p);

% G_eval = full(G_fun(x_k,p0))*0;
% H_eval = full(H_fun(x_k,p0))*0;
if settings.initial_comp_all_zero
    G_eval = zeros(n_comp,1);
    H_eval = G_eval;
else
    G_eval = full(G_fun(x_k,p0));
    H_eval = full(H_fun(x_k,p0));
end

% 0 \leq G(x) \perp H(x) \geq 0
ind_x1 = [];
ind_x2 = [];

% Lift complemetraties G
if settings.lift_complementarities_full
    % full lifting with duplicatse;
    x1 = SX.sym('x1',n_comp);
    lbx = [lbx;0*ones(n_comp,1)];
    ubx = [ubx;inf*ones(n_comp,1)];
    x = [x;x1];
    x_k = [x_k;G_eval];

    g = [g;x1-G];
    lbg = [lbg;0*ones(n_comp,1)];
    ubg = [ubg;0*ones(n_comp,1)];

    ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
    [ind_x1,~] = find(sparse(ind_x1_fun(x_k)==1));
elseif settings.lift_complementarities
    % lifting with only those that are not scaler
    [ind_scalar,ind_nonscalar, ind_map] = find_nonscalar(G,x);
    n_lift_x1 = length(ind_nonscalar);
    if n_lift_x1 == 0
        try x.jacobian(G_copy);
        catch
            n_lift_x1 = length(G_copy);
            ind_nonscalar = 1:n_lift_x1;
            ind_scalar = [];
        end
    end
    if n_lift_x1 > 0
        x1_lift = SX.sym('x1_lift',n_lift_x1);
        lbx = [lbx;0*ones(n_lift_x1,1)];
        ubx = [ubx;inf*ones(n_lift_x1 ,1)];
        x = [x;x1_lift];
        x_k = [x_k;G_eval(ind_nonscalar)];
        % x1 = [x(ind_scalar);x1_lift];
        % ind_x1 = [ind_map;n_primal_non_lifted+1:n_primal_non_lifted+n_lift_x1];
        % lift
        g = [g;x1_lift-G(ind_nonscalar)];
        lbg = [lbg;0*ones(n_lift_x1 ,1)];
        ubg = [ubg;0*ones(n_lift_x1 ,1)];

        x1 = G_copy;
        x1(ind_nonscalar) = x1_lift;
    else
        x1 = G;
        % ind_x1 = ind_map;
    end
else
    x1 = G; % (an expression, a nontivial function of x)
    if settings.update_comp_lower_bounds
        g = [g;G];
        lbg = [lbg;0*ones(n_comp,1)];
        ubg = [ubg; inf*ones(n_comp,1)];
    end
end

% Lift complemetraties H
if settings.lift_complementarities_full
    % full lifting with duplicatse;
    x2 = SX.sym('x2',n_comp);
    lbx = [lbx;0*ones(n_comp,1)];
    ubx = [ubx;inf*ones(n_comp,1)];
    x = [x;x2];
    x_k = [x_k;H_eval ];

    g = [g;x2-H];
    lbg = [lbg;0*ones(n_comp,1)];
    ubg = [ubg;0*ones(n_comp,1)];

    ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
    [ind_x2,~] = find(sparse(ind_x2_fun(x_k))==1);
elseif settings.lift_complementarities
    % lifting with only those that are not scaler
    [ind_scalar,ind_nonscalar, ind_map] = find_nonscalar(H,x);
    n_lift_x2 = length(ind_nonscalar);

    if n_lift_x2 == 0
        try x.jacobian(H_copy);
        catch
            n_lift_x2 = length(H_copy);
            ind_nonscalar = 1:n_lift_x2;
            ind_scalar = [];
        end
    end
    if n_lift_x2 > 0
        x2_lift = SX.sym('x2_lift',n_lift_x2);
        lbx = [lbx;0*ones(n_lift_x2,1)];
        ubx = [ubx;inf*ones(n_lift_x2 ,1)];
        x = [x;x2_lift];
        x_k = [x_k;H_eval(ind_nonscalar)];
        % ind_x2 = [ind_map;n_primal_non_lifted+n_lift_x1+1:n_primal_non_lifted+n_lift_x1+n_lift_x2];
        % lift
        g = [g;x2_lift-H(ind_nonscalar)];
        lbg = [lbg;0*ones(n_lift_x2 ,1)];
        ubg = [ubg;0*ones(n_lift_x2 ,1)];

        x2 = H_copy;
        x2(ind_nonscalar) = x2_lift;
    else
        x2 = H;
        % ind_x2 = ind_map;
    end
else
    x2 = H; % (an expression, a nontivial function of x)
    if settings.update_comp_lower_bounds
        g = [g;H];
        lbg = [lbg;0*ones(n_comp,1)];
        ubg = [ubg; inf*ones(n_comp,1)];
    end
end

if settings.lift_complementarities_full || settings.lift_complementarities
    ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
    [ind_x1,~] = find(sparse(ind_x1_fun(x_k)==1));
    ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
    [ind_x2,~] = find(sparse(ind_x2_fun(x_k))==1);
end
% update lb on x1 and x2
if settings.update_comp_lower_bounds
    lbx(ind_x1) = 0;
    lbx(ind_x2) = 0;
end
n_primal = length(x);
ind_x0 = [1:n_primal]';
ind_x0([ind_x1,ind_x2]) = [];


%% primal infeasiblity (% todo; there is repetition of code in creating the lpec)
% Split into equalites and inequalities
ind_g_eq = find(lbg == ubg);
ind_g_ineq = find(lbg < ubg);

ind_g_ineq_lb = find(lbg >- inf & lbg < ubg);
ind_g_ineq_ub = find(ubg < inf & lbg < ubg);

ind_x_lb = find(lbx > -inf);
ind_x_ub =  find(ubx < inf);

n_eq = length(ind_g_eq);
n_g_ineq_ub = length(ind_g_ineq_ub);
n_g_ineq_lb = length(ind_g_ineq_lb);

n_ubx = length(ind_x_ub);
n_lbx = length(ind_x_lb);

lbx_reduced = lbx(ind_x_lb);
ubx_reduced = ubx(ind_x_ub);
g_sym = g;
% Generate casadi functions for objective and constraint function evaluations
nabla_f = f.jacobian(x)';
% Zero order
g_eq = g(ind_g_eq)-lbg(ind_g_eq);                                          % g_eq = g - g_lb = 0
g_ineq_ub = ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
g_ineq_lb = g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0
g_ineq = [ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
    g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
n_ineq = size(g_ineq,1);

h_eq = max(abs(g_eq));
h_ineq_ub = max(min(g_ineq_ub,0));
h_ineq_lb = max(min(g_ineq_lb,0));
h_ubx = max(min(ubx-x,0));
h_lbx = max(min(x-lbx,0));
% Summary
h_std = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);


%% LPEC prepration
if settings.check_B_stationarity
    % Split into equalites and inequalities
    ind_g_eq = find(lbg == ubg);
    ind_g_ineq = find(lbg < ubg);

    ind_g_ineq_lb = find(lbg >- inf & lbg < ubg);
    ind_g_ineq_ub = find(ubg < inf & lbg < ubg);

    ind_x_lb = find(lbx > -inf);
    ind_x_ub =  find(ubx < inf);

    n_eq = length(ind_g_eq);
    n_g_ineq_ub = length(ind_g_ineq_ub);
    n_g_ineq_lb = length(ind_g_ineq_lb);

    n_ubx = length(ind_x_ub);
    n_lbx = length(ind_x_lb);

    lbx_reduced = lbx(ind_x_lb);
    ubx_reduced = ubx(ind_x_ub);
    g_sym = g;
    % Generate casadi functions for objective and constraint function evaluations
    nabla_f = f.jacobian(x)';
    % Zero order
    g_eq = g(ind_g_eq)-lbg(ind_g_eq);                                          % g_eq = g - g_lb = 0
    g_ineq_ub = ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
    g_ineq_lb = g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0
    g_ineq = [ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
        g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
    n_ineq = size(g_ineq,1);
    % first-order Constraint Jacobians
    if n_eq > 0
        nabla_g = g_sym.jacobian(x);
        nabla_g_eq = g_eq.jacobian(x);
    else
        nabla_g = [];
        nabla_g_eq = [];
    end
    if n_ineq > 0
        nabla_g_ineq_ub = g_ineq_ub.jacobian(x);
        nabla_g_ineq_lb = g_ineq_lb.jacobian(x);
    else
        nabla_g_ineq_ub = [];
        nabla_g_ineq_lb = [];
    end
    nabla_g_ineq = [nabla_g_ineq_ub; nabla_g_ineq_lb];

    % CasADi functions of function evaluations and derivaties
    % Zero order (objective and constraint function evaluations)
    f_fun =  Function('f_fun',{x,p},{f});
    g_eq_fun =  Function('g_eq_fun',{x,p},{g_eq});
    g_ineq_ub_fun = Function('g_ineq_ub_fun',{x,p},{g_ineq_ub});
    g_ineq_lb_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq_lb});
    g_ineq_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq});
    % First order (Gradients and Jacobian)
    nabla_f_fun = Function('nabla_f_fun',{x,p},{nabla_f});
    nabla_g_fun = Function('nabla_g_fun',{x,p},{nabla_g});
    nabla_g_eq_fun = Function('nabla_g_eq_fun',{x,p},{nabla_g_eq});
    nabla_g_ineq_ub_fun  = Function('nabla_g_ineq_ub_fun',{x,p},{nabla_g_ineq_ub});
    nabla_g_ineq_lb_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq_lb});
    nabla_g_ineq_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq});

    dims.n_primal = n_primal;
    dims.n_comp = n_comp;
    dims.n_eq = n_eq;
    dims.n_ineq = n_ineq;

    dims.ind_x0 = ind_x0;
    dims.ind_x1 = ind_x1;
    dims.ind_x2 = ind_x2;
    dims.n_slacks = 0;

    M = SX.sym('M', 1);
    y = SX.sym('y', dims.n_comp); % binary variablkes for comp. constraints
    % Big M reformulation of complementarities
    A_lpec_sym = [-x1+M*y; -x2-M*y];
    A_lpec = A_lpec_sym.jacobian([x;y]);
    A_lpec_fun = Function('A_lpec_fun',{M},{A_lpec});
    b_res = [x1;x2-M];
    b_res_fun = Function('b_res_fun',{x,p,M},{b_res});
    % for initalzing of an lpec
    sense_B = repmat('>',1,2*n_comp);
    sense = [repmat('=',1,dims.n_eq), repmat('>',1,dims.n_ineq), sense_B];
    vtype = [repmat('C',1,n_primal), repmat('B',1,n_comp)];
    vtype_num = [repmat(0,1,n_primal), repmat(1,1,n_comp)];
    A_lpec_k = A_lpec;
    % b_lpec_k = full(b_res_fun(x_k,p0));
    dims.n_auxiliary = n_comp;
    lb_binary = 0*ones(n_comp,1);
    ub_binary = 1*ones(n_comp,1);
    lpec_functions.A_lpec_fun = A_lpec_fun;
    lpec_functions.b_res_fun = b_res_fun;
    lpec_functions.vtype = vtype;
    lpec_functions.vtype_num = vtype_num;
    lpec_functions.sense = sense;
    lpec_functions.lb_binary = lb_binary;
    lpec_functions.ub_binary = ub_binary;
    lpec_functions.g_eq_fun = g_eq_fun;
    lpec_functions.f_fun = f_fun;
    lpec_functions.g_ineq_fun = g_ineq_fun;
    lpec_functions.nabla_g_eq_fun = nabla_g_eq_fun;
    lpec_functions.nabla_g_ineq_fun = nabla_g_ineq_fun;
    lpec_functions.nabla_f_fun = nabla_f_fun;
    lpec_functions.lbx = lbx;
    lpec_functions.ubx = ubx;
end
% lpec_casadi = create_lpec_functions(mpec_casadi,dims,settings,solver_initalization);


%% prepare relaxation method
switch settings.homotopy_parameter_steering
    case 'Direct_eq'
        if settings.aggregate_comps
            g_comp = x1'*x2-sigma*n_comp;
            lbg_comp = 0*ones(1,1);
            ubg_comp = 0*ones(1,1);
        else
            g_comp = x1.*x2-sigma;
            lbg_comp = 0*ones(n_comp,1);
            ubg_comp = 0*ones(n_comp,1);
        end
    case 'Direct'
        if settings.aggregate_comps
            g_comp = x1'*x2-sigma*n_comp;
            lbg_comp = -inf*ones(1,1);
            ubg_comp = 0*ones(1,1);
        else
            g_comp = x1.*x2-sigma;
            lbg_comp = -inf*ones(n_comp,1);
            ubg_comp = 0*ones(n_comp,1);
        end
    case 'Ell_1'
        g_comp = [];
        lbg_comp = [];
        ubg_comp = [];
        f = f+(x1'*x2)*(sigma)^(-1);
    case 'Ell_inf'
        s_eleastic = SX.sym('s_eleastic',1);
        f = f+(s_eleastic)*(sigma)^(-1);
        x = [x;s_eleastic];
        lbx = [lbx;0];
        ubx = [ubx;max(10,settings.sigma0*10)];
        x_k = [x_k;settings.sigma0];

        if settings.aggregate_comps
            g_comp = x1'*x2-s_eleastic*n_comp;
            lbg_comp = -inf*ones(1,1);
            ubg_comp = 0*ones(1,1);
        else
            g_comp = x1.*x2-s_eleastic;
            lbg_comp = -inf*ones(n_comp,1);
            ubg_comp = 0*ones(n_comp,1);
        end
    case 'None'
        settings.max_iter = 1;
        if settings.aggregate_comps
            g_comp = x1'*x2;
            lbg_comp = -inf*ones(1,1);
            ubg_comp = 0*ones(1,1);
        else
            g_comp = x1.*x2;
            lbg_comp = -inf*ones(n_comp,1);
            ubg_comp = 0*ones(n_comp,1);
        end
end

f_fun =  Function('f_fun',{x,p},{f});
h_std_fun  = Function('h_std_fun',{x,p},{h_std});

ind_comp = length(g)+1;
g = [g;g_comp];
lbg = [lbg; lbg_comp];
ubg = [ubg; ubg_comp];
ind_comp = ind_comp:1:length(g);

h_comp_constraints_tol_fun = Function('h_comp_constraints_tol_fun',{x,p},{max(min(abs(x1),abs(x2)))});
if settings.comp_res_bilinear
    % h_comp_constraints_fun = Function('h_comp_constraints_fun',{x,p},{max(abs(x1.*x2))});
    h_comp_constraints_fun = Function('h_comp_constraints_fun',{x,p},{max(abs(x1).*abs(x2))});
else
    h_comp_constraints_fun = Function('h_comp_constraints_fun',{x,p},{max(min(abs(x1),abs(x2)))});
end

%% nlp solver
nlp = struct('x', x, 'f', f, 'g', g,'p',p);
solver = nlpsol('solver', 'ipopt', nlp, settings.settings_casadi_nlp);
n_nlp_total = 0;
%% main homotopy loop
cpu_time_iters  = [];
success  = false;
problem_infeasible = false;
max_iterations_reached = false;

if settings.verbose_solver
    fprintf('iter\t obj\t comp_res\t inf_pr\t\t inf_du\t\t nlp_iters\t\t tau \t\t status\n')
    fprintf('---------------------------------------------------------------------------------------------------\n')
end

X_outer = [x_k];
cpu_time_nlp_iter = [];
%% Main homotopy loop
t_phase_ii_star = tic;

comp_res_ii =  full(h_comp_constraints_fun(x_k,p0));
f_k = full(f_fun(x_k,p0));
inf_pr_ii = full(h_std_fun(x_k,p0));
if settings.verbose_solver
    fprintf('%d \t %2.2e\t %2.2e\t  %2.2e\t  %2.2e\t %d \t\t %2.2e\t \t\t  %s  \n',0,f_k,comp_res_ii,inf_pr_ii,0,0,p0(end),'Inital guess');
end
comp_res_ii = 1e3;% avoid termination at feasible but not optimal point
while (comp_res_ii >= settings.comp_tol && ii <= settings.max_iter)
    t_nlp_start_ii = tic;
    solution = solver('x0',x_k,'p',p0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
    t_nlp_end_ii = toc(t_nlp_start_ii);
    n_nlp_total  = n_nlp_total + 1;
    cpu_time_nlp_iter = [cpu_time_nlp_iter,t_nlp_end_ii];
    stats = solver.stats();
    x_k = full(solution.x);
    comp_res_ii =  full(h_comp_constraints_fun(x_k,p0));
    f_k = full(f_fun(x_k,p0));
    inf_pr_ii = stats.iterations.inf_pr(end);
    inf_du_ii = stats.iterations.inf_du(end);
    if settings.verbose_solver
        fprintf('%d \t %2.2e\t %2.2e\t  %2.2e\t  %2.2e\t %d \t\t %2.2e\t \t\t  %s  \n',ii,f_k,comp_res_ii,inf_pr_ii,inf_du_ii,stats.iter_count,p0(end),(stats.return_status));
    end
    p0(end) = settings.kappa*p0(end);
    ii = ii+1;
    X_outer = [X_outer, x_k];
    if isequal(stats.return_status,'Infeasible_Problem_Detected') || isequal(stats.return_status,'Error_In_Step_Computation')
        break;
    end
end

%% status of last iter
if isequal(stats.return_status,'Solve_Succeeded') || isequal(stats.return_status,'Solved_To_Acceptable_Level') || isequal(stats.return_status,'Search_Direction_Becomes_Too_Small')
    success = true;
    solver_message = 'Stationary point was found successfuly.';
else
    success = false;
    solver_message = stats.return_status;
end

if isequal(stats.return_status,'Maximum_Iterations_Exceeded')
    success = false;
    solver_message = 'Last NLP in homotopy reached maximum number of iterations.';
end

if comp_res_ii > settings.comp_tol
    if ii == settings.max_iter
        max_iterations_reached = true;
        solver_message = 'Maximum number of iterations reached.';
    end
    success = false;
    solver_message = 'Complementarity tolerance not satisfied.';
end

if isequal(stats.return_status,'Infeasible_Problem_Detected')
    problem_infeasible  = true;
    success = false;
    solver_message = 'Infeasible problem detected.';
elseif isequal(stats.return_status,'Restoration_Failed')
    problem_infeasible  = true;
    success = false;
    solver_message = 'Infeasible problem detected.';
else
    problem_infeasible  = false;
end
cpu_time_phase_ii = toc(t_phase_ii_star);

%% biactive
G_res = full(G_fun(x_k(1:n_primal_non_lifted),p0));
H_res = full(H_fun(x_k(1:n_primal_non_lifted),p0));

if settings.check_B_stationarity
    active_set_estimate_k = find_active_sets(x_k, dims, settings.tol_active);
    n_biactive = sum(active_set_estimate_k.I_00);
else
    n_biactive = sum(G_res+H_res<settings.tol_active);
end

solution.G = G_res;
solution.H = H_res;
%% determine type of stationary point

multiplier_based_stationarity = 'X';
b_stationarity = false;
dims.ind_x1 = ind_x1;
dims.ind_x2 = ind_x2;
dims.n_primal = n_primal;
dims.n_comp = n_comp;

x_k_before_tnlp = x_k;
%%
% ----------------------------------- MULTIPLIER BASED STAT -------------------------------------
N_TNLP = 3; % max attemps to solve a TNLP;
ii = 1;
tol_ell_inf = full(h_comp_constraints_tol_fun(x_k,p0));
% tol_active = max(2*settings.tol_active*1e3,2*tol_ell_inf);
tol_active = settings.tol_active;
if settings.compute_tnlp_stationary_point && success && settings.lift_complementarities && ~settings.problem_is_lpec
    fprintf('----------------------------------- determining stationarity -----------------------------------------------\n')

    % Solve a TNLP:
    ubg(ind_comp) = +inf;  % remove upper bound on bilinear terms to avodi multipliers
    if strcmp(settings.homotopy_parameter_steering,'Ell_1') || strcmp(settings.homotopy_parameter_steering,'Ell_inf')
        ubx(end) = +inf;
    end
    while  ii <= N_TNLP
        active_set_estimate_k = find_active_sets(x_k, dims, tol_active);
        ubx(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
        ubx(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
        ubx(dims.ind_x1(active_set_estimate_k.I_00)) = 0;
        ubx(dims.ind_x2(active_set_estimate_k.I_00)) = 0;
        % n_biactive = sum(active_set_estimate_k.I_00)

        solution = solver('x0',x_k,'p',p0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
        lambda_x_tnlp = full(solution.lam_x);
        x_tnlp = full(solution.x);
        f_tnlp = full(solution.f);
        comp_res_tnlp =  full(h_comp_constraints_fun(x_k,p0));
        inf_pr_tnlp = stats.iterations.inf_pr(end);
        inf_du_tnlp = stats.iterations.inf_du(end);
        settings.tol_active = tol_active;
        if settings.verbose_solver
            fprintf('%d \t %2.2e\t %2.2e\t  %2.2e\t  %2.2e\t %d \t\t %2.2e\t \t\t  %s  \n',ii,f_tnlp,comp_res_tnlp,inf_pr_tnlp,inf_du_tnlp,stats.iter_count,0,(stats.return_status));
        end
        % fprintf('\t\t ||x_tnlp - x_k|| = %2.2e, |f_tnlp-f_k| = %2.2e \n', norm(x_tnlp-x_k,inf),abs(f_tnlp-f_k))
        if norm(x_tnlp-x_k,inf) <= 1e-6 || abs(f_tnlp-f_k)/abs(f_k) <= 1e-3 || ii == N_TNLP
            [multiplier_based_stationarity, ~] = determine_multipliers_based_stationary_point(x_tnlp,lambda_x_tnlp,dims,settings);
            n_biactive = sum(active_set_estimate_k.I_00);
            if ii~=N_TNLP
                x_k = x_tnlp;
            end
            break;
        else
            tol_active = 0.1*tol_active;
            ii = ii+1;
        end
    end
    if ~strcmp(multiplier_based_stationarity,'X')
    end
end

% if strcmp(multiplier_based_stationarity,'W')
%     keyboard;
% end
%% Check B stationarity
% ----------------------------------- B STAT -------------------------------------
f_lpec = 1;
rho_TR = 1e-3;

if ~settings.problem_is_lpec && success
    settings_lpec = LPECSolverOptions();
    N_LPEC = 3;
    ii = 1;
    lpec = create_lpec_subproblem(x_k(1:n_primal),p0 ,rho_TR, lpec_functions, dims, settings, settings.tol_active);
    nabla_f_k = full(nabla_f_fun(x_k(1:n_primal),p0));
    if norm(nabla_f_k) >= 1e2 && settings.rescale_large_objective_gradients
        nabla_f_k = nabla_f_k./(norm(nabla_f_k));
        lpec.f = nabla_f_k;
    end
    while ii <= N_LPEC
        % full(h_comp_constraints_fun(x_k,p0))
        lpec.rho_TR = rho_TR;
        [results_lpec,stats_lpec] = lpec_solver(lpec,settings_lpec);
        f_lpec = results_lpec.f_opt;
        d_lpec = results_lpec.d_lpec;
        if abs(f_lpec) <= settings.tol_B_stationarity
            % if norm((d_lpec),inf)<= settings.tol_B_stationarity
            b_stationarity = true;
            break;
        end
        rho_TR = 0.1*rho_TR;
        ii = ii+1;
    end
    f_lpec = results_lpec.f_opt;
end

% double check biactive set
if ~settings.problem_is_lpec && b_stationarity && success  && ~strcmp(multiplier_based_stationarity,'S')
    % active_set_estimate_k = find_active_sets(x_k(1:n_primal)+results_lpec.d_lpec, dims, settings.tol_active);
    active_set_estimate_k = find_active_sets(x_tnlp(1:n_primal), dims, settings.tol_active);
    n_biactive = sum(active_set_estimate_k.I_00);
    if n_biactive == 0
        multiplier_based_stationarity = 'S';
    end
end

% use LPEC active set, if difference not too big take this as solution
if ~settings.problem_is_lpec && success && ~b_stationarity
    % if strcmp(multiplier_based_stationarity,'X')
    try
        % keyboard;
        active_set_estimate_k = find_active_sets(x_k+results_lpec.d_lpec, dims, settings.tol_active);
        ubx(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
        ubx(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
        ubx(dims.ind_x1(active_set_estimate_k.I_00)) = 0;
        ubx(dims.ind_x2(active_set_estimate_k.I_00)) = 0;
        
        solution = solver('x0',x_k,'p',p0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
        x_tnlp = full(solution.x);
        f_tnlp = full(solution.f);
        lambda_x_tnlp = full(solution.lam_x);
        comp_res_tnlp =  full(h_comp_constraints_fun(x_k,p0));
        inf_pr_tnlp = stats.iterations.inf_pr(end);
        inf_du_tnlp = stats.iterations.inf_du(end);
        fprintf('----------------------------------- lpec active set -----------------------------------------------\n')
        if settings.verbose_solver
            fprintf('%d \t %2.2e\t %2.2e\t  %2.2e\t  %2.2e\t %d \t\t %2.2e\t \t\t  %s  \n',ii,f_tnlp,comp_res_tnlp,inf_pr_tnlp,inf_du_tnlp,stats.iter_count,0,(stats.return_status));
        end
        fprintf('\t\t ||x_tnlp - x_k|| = %2.2e, |f_tnlp-f_k| = %2.2e \n', norm(x_tnlp-x_k,inf),abs(f_tnlp-f_k))
        % Compute multipliers
        ii = 1;
        rho_TR = 1e-3;
        if norm(x_tnlp-x_k,inf) <= 1e-6 || abs(f_tnlp-f_k)/abs(f_k+1e-16) <= 1e-6
            [multiplier_based_stationarity, ~] = determine_multipliers_based_stationary_point(x_tnlp,lambda_x_tnlp,dims,settings);
            active_set_estimate_k = find_active_sets(x_tnlp, dims, settings.tol_active);
            n_biactive = sum(active_set_estimate_k.I_00);
            % check b stationarity
            lpec = create_lpec_subproblem(x_tnlp(1:n_primal),p0 ,rho_TR, lpec_functions, dims, settings, settings.tol_active);
        else
            lpec = create_lpec_subproblem(x_k_before_tnlp(1:n_primal),p0 ,rho_TR, lpec_functions, dims, settings, settings.tol_active);
        end
        while ii <= N_LPEC
            lpec.rho_TR = rho_TR;
            [results_lpec,stats_lpec] = lpec_solver(lpec,settings_lpec);
            f_lpec = results_lpec.f_opt;
            d_lpec = results_lpec.d_lpec;
            % if abs(f_lpec) <= settings.tol_B_stationarity
            if abs(f_lpec) <= settings.tol_B_stationarity
                b_stationarity = true;
                x_k = x_tnlp;
                break;
            end
            rho_TR = 0.1*rho_TR;
            ii = ii+1;
        end
    catch
    end
end

if strcmp(multiplier_based_stationarity,'S')
    b_stationarity = true;
end
if settings.verbose_summary
    fprintf('\n');
    print_iter_summary_scholtes(f_k,inf_pr_ii,comp_res_ii,solver_message,multiplier_based_stationarity,b_stationarity,n_biactive,f_lpec,rho_TR);
end

% if strcmp(multiplier_based_stationarity,'X')
%     keyboard;
% end

% if ~b_stationarity
%     keyboard;
% end
% if ~b_stationarity && success && ~settings.problem_is_lpec
%     f_lpec
%     keyboard;
% end
%% results
solution.f = f_k;
solution.x_lifted = x_k;
solution.x = x_k(1:n_primal_non_lifted);
solution.x1  = full(G_fun(x_k(1:n_primal_non_lifted),p0));
solution.x2 =  full(H_fun(x_k(1:n_primal_non_lifted),p0));

stats.comp_res = comp_res_ii;
stats.success = success;
stats.success_phase_i = true;
stats.f_lpec = f_lpec;
stats.problem_infeasible  = problem_infeasible;
stats.max_iterations_reached = max_iterations_reached;
stats.solved_in_phase_i = false; % homotopy sovlers dont have a phase i
stats.multiplier_based_stationarity = multiplier_based_stationarity;
stats.b_stationarity = b_stationarity;
stats.n_biactive = n_biactive;
stats.success_despite_infeasiblity = success;
% total time per phase
stats.cpu_time_total = cpu_time_phase_ii;
stats.cpu_time_phase_i = 0;
stats.cpu_time_phase_ii = cpu_time_phase_ii;
% nlp time per phase
stats.cpu_time_nlp = sum(cpu_time_nlp_iter);
stats.cpu_time_nlp_iter = cpu_time_nlp_iter;
stats.cpu_time_nlp_phase_i = 0;
stats.cpu_time_nlp_phase_ii = sum(cpu_time_nlp_iter);
stats.n_nlp_total = n_nlp_total;
%lpec time
stats.cpu_time_lpec = 0;
stats.cpu_time_lpec_phase_i = 0;
stats.cpu_time_lpec_phase_ii = 0;

% cpu time per iter
stats.iter.cpu_time_lpec_phase_i_iter = 0;
stats.iter.cpu_time_lpec_phase_ii_iter = 0;
stats.iter.cpu_time_nlp_phase_i_iter = 0;
stats.iter.cpu_time_nlp_phase_ii_iter = cpu_time_nlp_iter;

% Dummy Lpec data
stats.n_lpec_total = 0;
stats.iter.nodecount_phase_i = 0;
stats.iter.nodecount_phase_ii= 0;
stats.iter.baritercount_phase_i = 0;
stats.iter.baritercount_phase_ii= 0;
stats.iter.itercount_phase_i = 0;
stats.iter.itercount_phase_ii= 0;

if sum(cpu_time_nlp_iter)>0
    stats.iter.cpu_time_nlp_iter = cpu_time_nlp_iter;
else
    stats.iter.cpu_time_nlp_iter = nan+cpu_time_nlp_iter;
end
stats.iter.active_set_changes = 0;
stats.iter.X_outer = X_outer;
end

