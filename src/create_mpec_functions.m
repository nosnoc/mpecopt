function [mpec_casadi, dims, solver_initalization, stats] =  create_mpec_functions(mpec,solver_initalization,settings)
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

x_k = solver_initalization.x0; % inital guess

if isfield(mpec,'p')
    p = mpec.p;
    if ~isfield(solver_initalization,'p0')
        error('mpec.p defined, but no solver_initalization.p0 passed.')
    end
else
    p = [];
    solver_initalization.p0 = [];
end

%% Edit complementarity constraints
n_primal_non_lifted = length(solver_initalization.x0);
n_comp = size(G,1);
G_fun = Function('G_fun',{x,p},{G});
H_fun = Function('H_fun',{x,p},{H});

G_copy = G_fun(x,p);
H_copy = H_fun(x,p);


if settings.initial_comp_all_zero
    G_eval = zeros(n_comp,1);
    H_eval = zeros(n_comp,1);
else
    G_eval = full(G_fun(x_k,solver_initalization.p0));
    H_eval = full(G_fun(x_k,solver_initalization.p0));
end




% check if the comps are expressions or just subvectors of w.
if settings.lift_complementarities_full
    % define lift vairables
    x1 = SX.sym('x1',n_comp);
    solver_initalization.lbx = [solver_initalization.lbx;0*ones(n_comp,1)];
    solver_initalization.ubx = [solver_initalization.ubx;inf*ones(n_comp,1)];
    % update x and init guess
    x = [x;x1];
    x_k = [x_k;G_eval];
    % lift
    g = [g;x1-G];
    solver_initalization.lbg = [solver_initalization.lbg;0*ones(n_comp,1)];
    solver_initalization.ubg = [solver_initalization.ubg;0*ones(n_comp,1)];
else
    % lifting with only those that are not scaler
    % define lift vairables
    [ind_scalar,ind_nonscalar_x1, ind_map] = find_nonscalar(G,x);
    n_lift_x1 = length(ind_nonscalar_x1);
    if n_lift_x1 == 0
        try x.jacobian(G_copy);
        catch
            n_lift_x1 = length(G_copy);
            ind_nonscalar_x1 = 1:n_lift_x1;
            ind_scalar = [];
        end
    end
    if n_lift_x1 > 0
        x1_lift = SX.sym('x1_lift',n_lift_x1);
        solver_initalization.lbx = [solver_initalization.lbx;0*ones(n_lift_x1,1)];
        solver_initalization.ubx = [solver_initalization.ubx;inf*ones(n_lift_x1 ,1)];
        % x1 = [x(ind_scalar);x1_lift];
        % update x and init guess
        x = [x;x1_lift];
        x_k = [x_k;G_eval(ind_nonscalar_x1)];
        % lift
        g = [g;x1_lift-G(ind_nonscalar_x1)];
        solver_initalization.lbg = [solver_initalization.lbg;0*ones(n_lift_x1 ,1)];
        solver_initalization.ubg = [solver_initalization.ubg;0*ones(n_lift_x1 ,1)];

        x1 = G_copy;
        x1(ind_nonscalar_x1) = x1_lift;
    else
        x1 = G;
    end
end

if settings.lift_complementarities_full
    % define lift vairables
    x2 = SX.sym('x2',n_comp);
    solver_initalization.lbx = [solver_initalization.lbx;0*ones(n_comp,1)];
    solver_initalization.ubx = [solver_initalization.ubx;inf*ones(n_comp,1)];
    % update x and init guess
    x = [x;x2];
    x_k = [x_k;H_eval];
    % lift
    g = [g;x2-H];
    solver_initalization.lbg = [solver_initalization.lbg;0*ones(n_comp,1)];
    solver_initalization.ubg = [solver_initalization.ubg;0*ones(n_comp,1)];
else
    % lifting with only those that are not scaler
    [ind_scalar,ind_nonscalar_x2, ind_map] = find_nonscalar(H,x);
    n_lift_x2 = length(ind_nonscalar_x2);
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
        solver_initalization.lbx = [solver_initalization.lbx;0*ones(n_lift_x2,1)];
        solver_initalization.ubx = [solver_initalization.ubx;inf*ones(n_lift_x2 ,1)];
        % x2 = [x(ind_scalar);x2_lift];
        % update x and init guess
        x = [x;x2_lift];
        x_k = [x_k;H_eval(ind_nonscalar_x2)];
        % lift
        g = [g;x2_lift-H(ind_nonscalar_x2)];
        solver_initalization.lbg = [solver_initalization.lbg;0*ones(n_lift_x2 ,1)];
        solver_initalization.ubg = [solver_initalization.ubg;0*ones(n_lift_x2 ,1)];

        x2 = H_copy;
        x2(ind_nonscalar_x2) = x2_lift;
    else
        x2 = H;
    end
end

% find index set
if n_comp > 0
    ind_x1 = [];
    ind_x2 = [];
    ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
    [ind_x1,~] = find(sparse(ind_x1_fun(x_k)==1));
    ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
    [ind_x2,~] = find(sparse(ind_x2_fun(x_k))==1);
    settings.nlp_is_mpec = 1; % TODO: check is this still used? (its purpose: if not an mpec, just make single nlp call without mpec machinery);
else
    settings.nlp_is_mpec = 0;
    ind_x1 = [];
    ind_x2 = [];
end

solver_initalization.lbx(ind_x1) = 0;
solver_initalization.lbx(ind_x2) = 0;

n_primal = length(x);
n_primal_x0 = n_primal - 2*n_comp; % primal variables excluding the complementarity variables;
ind_x0 = [1:n_primal]';
if settings.nlp_is_mpec
    ind_x0([ind_x1,ind_x2]) = [];
end
x0 = x(ind_x0); % Variables not involved in complementarity constraints.

%% Split into equalites and inequalities
ind_g_eq = find(solver_initalization.lbg == solver_initalization.ubg);
ind_g_ineq = find(solver_initalization.lbg < solver_initalization.ubg);

ind_g_ineq_lb = find(solver_initalization.lbg >- inf & solver_initalization.lbg < solver_initalization.ubg);
ind_g_ineq_ub = find(solver_initalization.ubg < inf & solver_initalization.lbg < solver_initalization.ubg);

ind_x_lb = find(solver_initalization.lbx > -inf);
ind_x_ub =  find(solver_initalization.ubx < inf);

n_eq = length(ind_g_eq);
n_g_ineq_ub = length(ind_g_ineq_ub);
n_g_ineq_lb = length(ind_g_ineq_lb);

n_ubx = length(ind_x_ub);
n_lbx = length(ind_x_lb);

g_sym = g;
% Generate casadi functions for objective and constraint function evaluations
nabla_f = f.jacobian(x)';
% Zero order
g_eq = g(ind_g_eq)-solver_initalization.lbg(ind_g_eq);                                          % g_eq = g - g_lb = 0
g_ineq_ub = solver_initalization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
g_ineq_lb = g(ind_g_ineq_lb)-solver_initalization.lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0
g_ineq = [solver_initalization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
    g(ind_g_ineq_lb)-solver_initalization.lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
n_ineq = size(g_ineq,1);

% g_tnlp = [g_eq;g_ineq];
% lbg_tnlp = [zeros(n_eq,1);zeros(n_ineq,1)];
% ubg_tnlp = [zeros(n_eq,1);inf*ones;zeros(n_ineq,1)];

% first-order constraint Jacobians
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

%% CasADi functions for constraint evaluations and their derivaties
% Zero order (objective and constraint function evaluations)
mpec_casadi = struct;

mpec_casadi.f = f;
mpec_casadi.g = g;
mpec_casadi.x = x;
mpec_casadi.x1 = x1;
mpec_casadi.x2 = x2;
mpec_casadi.x0 = x0;
mpec_casadi.p = p;

mpec_casadi.f_fun =  Function('f_fun',{x,p},{f});
mpec_casadi.g_eq_fun =  Function('g_eq_fun',{x,p},{g_eq});
mpec_casadi.g_ineq_ub_fun = Function('g_ineq_ub_fun',{x,p},{g_ineq_ub});
mpec_casadi.g_ineq_lb_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq_lb});
mpec_casadi.g_ineq_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq});
mpec_casadi.G_fun = G_fun;
mpec_casadi.H_fun = H_fun;
% First order (Gradients and Jacobian)
mpec_casadi.nabla_f_fun = Function('nabla_f_fun',{x,p},{nabla_f});
mpec_casadi.nabla_g_fun = Function('nabla_g_fun',{x,p},{nabla_g});
mpec_casadi.nabla_g_eq_fun = Function('nabla_g_eq_fun',{x,p},{nabla_g_eq});
mpec_casadi.nabla_g_ineq_ub_fun  = Function('nabla_g_ineq_ub_fun',{x,p},{nabla_g_ineq_ub});
mpec_casadi.nabla_g_ineq_lb_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq_lb});
mpec_casadi.nabla_g_ineq_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq});


%% Infeasiblity meausre in inf norm
% all standard constraints
h_eq = max(abs(g_eq));
h_ineq_ub = max(min(g_ineq_ub,0));
h_ineq_lb = max(min(g_ineq_lb,0));
h_ubx = max(min(solver_initalization.ubx-x,0));
h_lbx = max(min(x-solver_initalization.lbx,0));
% Summary
h_std = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);
if n_comp > 0
    if settings.comp_res_bilinear
        h_comp= max(abs(x1).*abs(x2)); % here kappa =0.1 to have reduce of the value by factor of 10
        % h_comp= sqrt(max(abs(min((x1),(x2))))); % here kappa =0.01 to reduce the value above by factor of 10
    else
        h_comp= max(min(abs(x1),abs(x2))); % here kappa =0.01 to reduce the value above by factor of 10
        % h_comp = max(abs(min(x1,x2)));
    end
else
    h_comp  = 0;
end

%% CasADi Functions for constraints infeasiblity
mpec_casadi.h_std_fun  = Function('h_std_fun',{x,p},{h_std});
mpec_casadi.h_comp_fun  = Function('h_comp_fun',{x,p},{h_comp});
mpec_casadi.h_total_fun = Function('h_comp_fun',{x,p},{max(h_comp,h_std)});

solver_initalization.x0 = x_k;
%% Store some dimensions
dims.n_slacks = 0; % in generla no slacks, except in feasiblity problems
dims.ind_x0 = ind_x0;
dims.ind_x1 = ind_x1;
dims.ind_x2 = ind_x2;
dims.n_primal = n_primal;
dims.n_comp = n_comp;
dims.n_eq = n_eq;
dims.n_ineq = n_ineq;
dims.n_primal_non_lifted = n_primal_non_lifted;
dims.n_lift_x1 = n_lift_x1;
dims.n_lift_x2 = n_lift_x2;
dims.n_auxiliary = dims.n_comp; % number of binary variables in LPEC

% index sets of general constraints
dims.ind_g_eq = ind_g_eq;
dims.ind_g_ineq = ind_g_ineq;
dims.ind_g_ineq_lb = ind_g_ineq_lb;
dims.ind_g_ineq_ub = ind_g_ineq_lb; % if the last two have same indicies then it is a two sided ineq;
dims.ind_nonscalar_x1 = ind_nonscalar_x1;
dims.ind_nonscalar_x2 = ind_nonscalar_x2;

%%  Main piece NLP solver (BNLP or TNLP)
t_generate_nlp_solvers = tic;
nlp = struct('x', x,'f', f,'g', g,'p',p);
solver = nlpsol('solver', 'ipopt', nlp, settings.settings_casadi_nlp);
stats.cpu_time_generate_nlp_solvers = toc(t_generate_nlp_solvers);
mpec_casadi.solver = solver;

end

