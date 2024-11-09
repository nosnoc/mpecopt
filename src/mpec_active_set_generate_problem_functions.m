import casadi.*
%% Generate casadi functions for objective and constraint function evaluations
nabla_f = f.jacobian(x)';
% Zero order
g_eq = g(ind_g_eq)-lbg(ind_g_eq);                                          % g_eq = g - g_lb = 0
g_ineq_ub = ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
g_ineq_lb = g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0

g_ineq = [ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
    g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
n_ineq = size(g_ineq,1);

g_tnlp = [g_eq;g_ineq];
lbg_tnlp = [zeros(n_eq,1);zeros(n_ineq,1)];
ubg_tnlp = [zeros(n_eq,1);inf*ones;zeros(n_ineq,1)];

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

%% CasADi functions of function evaluations and derivaties
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

%% Infeasiblity meausre in inf norm
% all standard constraints
h_eq = max(abs(g_eq));
h_ineq_ub = max(min(g_ineq_ub,0));
h_ineq_lb = max(min(g_ineq_lb,0));
h_ubx = max(min(ubx-x,0));
h_lbx = max(min(x-lbx,0));
% Summary
h_standard_constraints = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);
if n_comp > 0
    h_comp_constraints = max(abs(min(x1,x2)));
else
    h_comp_constraints  = 0;
end

%% CasADi Functions for constraints infeasiblity
h_standard_constraints_fun  = Function('h_standard_constraints_fun',{x,p},{h_standard_constraints});
h_comp_constraints_fun  = Function('h_comp_constraints_fun',{x,p},{h_comp_constraints});
h_total_fun = Function('h_comp_constraints_fun',{x,p},{max(h_comp_constraints,h_standard_constraints)});

%% Store some dimensions
dims.ind_x0  = ind_x0;
dims.ind_x1  = ind_x1;
dims.ind_x2  = ind_x2;
dims.n_primal = n_primal;
dims.n_comp = n_comp;
dims.n_eq = n_eq;
dims.n_ineq = n_ineq;

