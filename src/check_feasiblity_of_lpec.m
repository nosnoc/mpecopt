function [lpec_is_feasible,min_rho_TR,lpec_infeasiblity] = check_feasiblity_of_lpec(lpec,settings)

x_lin = lpec.x_lin;
comp_res = max(min(abs(x_lin(lpec.dims.ind_x1)),abs(x_lin(lpec.dims.ind_x2))));
bounds_res = max([min(0,x_lin-lpec.lb);min(0,lpec.ub-x_lin)]);
eq_res = max(lpec.b_eq); % lpec.A_eq
ineq_res = min(0,lpec.b_ineq); 

lpec_infeasiblity = max([comp_res;bounds_res;eq_res;ineq_res]);

if lpec_infeasiblity > settings.tol
    lpec_is_feasible = false;
else
    lpec_is_feasible = true;
end

min_rho_TR = 10*comp_res+settings.tol_B_stationarity*2;


% A_gurobi = sparse([[lpec.A_eq, zeros(lpec.dims.n_eq, lpec.dims.n_auxiliary)];...
%             [lpec.A_ineq,zeros(lpec.dims.n_ineq, lpec.dims.n_auxiliary)];...
%             lpec.A_lpec]);
%         b_gurobi = [-lpec.b_eq; -lpec.b_ineq; lpec.b_lpec];
end