function [ocp_solver,stats] = mpec_active_set_nosnoc(ocp_solver,solver_settings)
% this function takes in a nosnoc discrete time problem and solves it with mpec active set
% extract mpec symbolic expressions
mpec = ocp_solver.discrete_time_problem.to_casadi_struct;
mpec.G= ocp_solver.discrete_time_problem.G.sym;
mpec.H= ocp_solver.discrete_time_problem.H.sym;
mpec.g= ocp_solver.discrete_time_problem.g.sym;
mpec.x= ocp_solver.discrete_time_problem.w.sym;
mpec.f= ocp_solver.discrete_time_problem.f;

% extract solver initalization
x0 = ocp_solver.discrete_time_problem.w.init;
lbx = ocp_solver.discrete_time_problem.w.lb;
ubx = ocp_solver.discrete_time_problem.w.ub;
lbg = ocp_solver.discrete_time_problem.g.lb;
ubg = ocp_solver.discrete_time_problem.g.ub;
p0 = ocp_solver.discrete_time_problem.p.val;

solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg,'p0',p0);
%% The problem
[solution,stats] = mpec_optimizer(mpec,solver_initalization,solver_settings);
w_opt = full(solution.x);
ocp_solver.discrete_time_problem.w.res  = w_opt;
end

