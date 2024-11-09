function  [lpec_casadi] = create_lpec_functions(mpec_casadi,dims,settings,solver_initalization)
import casadi.*

% symoblics
x1 = mpec_casadi.x1;
x2 = mpec_casadi.x2;
x = mpec_casadi.x;
% x1 = x(dims.ind_x1);
% x2 = x(dims.ind_x2);
p = mpec_casadi.p;
M = SX.sym('M', 1);
y = SX.sym('y', dims.n_comp); % binary variablkes for comp. constraints

% Big M reformulation of complementarities
A_lpec_sym = [-x1+M*y; -x2-M*y];
b_res = [x1;x2-M];
A_lpec = A_lpec_sym.jacobian([x;y]);
A_lpec_fun = Function('A_lpec_fun',{M},{A_lpec});
b_res_fun = Function('b_res_fun',{x,p,M},{b_res});

% for initalzing of an lpec
sense_B = repmat('>',1,2*dims.n_comp);
sense = [repmat('=',1,dims.n_eq), repmat('>',1,dims.n_ineq), sense_B];
vtype = [repmat('C',1,dims.n_primal), repmat('B',1,dims.n_comp)];
vtype_num = [zeros(1,dims.n_primal), ones(1,dims.n_comp)];

lb_binary = 0*ones(dims.n_comp,1);
ub_binary = 1*ones(dims.n_comp,1);

lpec_casadi.A_lpec_fun = A_lpec_fun;
lpec_casadi.b_res_fun = b_res_fun;
lpec_casadi.vtype = vtype;
lpec_casadi.vtype_num = vtype_num;
lpec_casadi.sense = sense;
lpec_casadi.lb_binary = lb_binary;
lpec_casadi.ub_binary = ub_binary;
lpec_casadi.lbx = solver_initalization.lbx;
lpec_casadi.ubx = solver_initalization.ubx;
% Copy mpec casadi functions needed for the linearization
lpec_casadi.g_eq_fun = mpec_casadi.g_eq_fun;
lpec_casadi.f_fun = mpec_casadi.f_fun;
lpec_casadi.g_ineq_fun = mpec_casadi.g_ineq_fun;
lpec_casadi.nabla_g_eq_fun = mpec_casadi.nabla_g_eq_fun;
lpec_casadi.nabla_g_ineq_fun = mpec_casadi.nabla_g_ineq_fun;
lpec_casadi.nabla_f_fun = mpec_casadi.nabla_f_fun;
end