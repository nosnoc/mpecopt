function x_k = project_to_bounds(x_k,lbx,ubx,dims)
    x_lin_0 = x_k(dims.ind_x0);%
    x_lin_1 = x_k(dims.ind_x1);%
    x_lin_2 = x_k(dims.ind_x2);%
    ind_x0_infeasible = x_lin_0<lbx(dims.ind_x0) | x_lin_0>ubx(dims.ind_x0);
    % if not in bounds put into middle, if bounds exist; if upper or lower bound -inf, set to zero;
    x_lin_0(ind_x0_infeasible) = 0.5*(max(lbx(ind_x0_infeasible),-1e3)+min(ubx(ind_x0_infeasible),1e3)); % mean between reduced bounds
    % project to satisfy nonnegativity
    x_lin_1 = max(0,x_lin_1);
    x_lin_2 = max(0,x_lin_2);
    % project to satisfy complementarity
    x_lin_1(x_lin_1<=x_lin_2) = 0;
    x_lin_2(x_lin_1>x_lin_2) = 0; % cf. Section 2.3 in Kirches2022
    x_k(dims.ind_x0) = x_lin_0;
    x_k(dims.ind_x1) = x_lin_1;
    x_k(dims.ind_x2) = x_lin_2;
end