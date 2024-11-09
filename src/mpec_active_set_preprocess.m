import casadi.*

%% check nlp input
if isfield(nlp,'x')
    x = vertcat(nlp.x{:});
    n_primal = length(x);
    if isfield(problem_data,'x0')
        x_k = problem_data.x0;
    else
        x_k = ones(n_primal,1);
    end
else
    error('Please provide the vector of degrees of freedom x.')
end

if isfield(nlp,'mpcc')
    mpcc = nlp.mpcc;
end

%% check nlp input
if isfield(nlp,'f')
    f = nlp.f;
else
    f = 0;
    warning('No objective was provided.');
end

if isfield(nlp,'g')
    g = nlp.g;
    n_constraints = length(g);
else
    n_constraints = 0;
end

if isfield(nlp,'comp1') && isfield(nlp,'comp2')
    x1 = nlp.comp1;
    x2 = nlp.comp2;
    n_comp = size(x1,1);
    if length(x1)~= length(x2)
        error('The vector comp1 and comp2 must have the same length.')
    end
    % find index set
    if n_comp > 0
        ind_x1 = [];
        ind_x2 = [];
        tic
        ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
        [ind_x1,~] = find(sparse(ind_x1_fun(x_k)==1));
        ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
        [ind_x2,~] = find(sparse(ind_x2_fun(x_k))==1);

        settings.nlp_is_mpec = 1;
    else
        settings.nlp_is_mpec = 0;
        ind_x1 = [];
        ind_x2 = [];
    end
else
    settings.nlp_is_mpec= 0;
end

n_primal_x0 = n_primal - 2*n_comp; % primal variables excluding the complementarity variables;
ind_x0 = [1:n_primal]';
if settings.nlp_is_mpec
    ind_x0([ind_x1,ind_x2]) = [];
end
% Variables not involved in complementarity constraints.
x0 = x(ind_x0);

%% Parameters
if isfield(nlp,'p')
    p = nlp.p;
    if isfield(problem_data,'p')
        p_val = problem_data.p;
        if length(p_val)~=length(p)
            error('Length of p and its value must be the same.')
        end
    end
else
    p = [];
    p_val = [];
end

%% check problem data input and read
if any(problem_data.lbg > problem_data.ubg)
    error('For some component it holds that lbg > ubg.');
else
    lbg = problem_data.lbg;
    ubg = problem_data.ubg;
end
if any(problem_data.lbx > problem_data.ubx)
    error('For some component it holds that lbx > ubx.');
else
    lbx = problem_data.lbx;
    ubx = problem_data.ubx;
end

if length(lbx)~= length(x) || length(ubx) ~= length(x)
    error('lbx, ubx and x must be vectors of the same length.')
end

if length(lbg)~= length(g) || length(ubg) ~= length(g)
    error('lbg, ubg and g must be vectors of the same length.')
end

% update complementarity lower bounds if they do not exist
if settings.nlp_is_mpec
    lbx(ind_x1 == -inf) = 0;
    lbx(ind_x2 == -inf) = 0;
end

%% Split into equalites and inequalities
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