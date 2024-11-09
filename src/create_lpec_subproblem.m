function lpec = create_lpec_subproblem(x_lin,p0,rho_TR,lpec_functions,dims,settings,tol_active_set)

% Default values
lpec.lb = lpec_functions.lbx;
lpec.ub = lpec_functions.ubx;
lpec.x_lin = x_lin; % linearization point
lpec.dims  = dims; %
lpec.sense = lpec_functions.sense;
lpec.vtype = lpec_functions.vtype;
lpec.vtype_num = lpec_functions.vtype_num;
lpec.lb_binary = lpec_functions.lb_binary; % relevant if the lpec is reduced
lpec.ub_binary = lpec_functions.ub_binary;

% Modifications

lpec.b_eq  = full(lpec_functions.g_eq_fun(x_lin,p0));
lpec.b_ineq  = full(lpec_functions.g_ineq_fun(x_lin,p0));
% f_k  = full(f_fun(x_lin,p0)); %% NEEDED outside of this function?
% Evaluate LPEC data - First-order
nabla_f  = full(lpec_functions.nabla_f_fun(x_lin,p0));
lpec.f = nabla_f;

lpec.A_eq  = full(lpec_functions.nabla_g_eq_fun(x_lin,p0));
lpec.A_ineq  = full(lpec_functions.nabla_g_ineq_fun(x_lin,p0));

M_val = max([settings.BigM, 1.01*max(x_lin([dims.ind_x1;dims.ind_x2])),rho_TR]);

if settings.consider_all_complementarities_in_lpec
    A_lpec = full(lpec_functions.A_lpec_fun(M_val));
    lpec.A_lpec = A_lpec;
    lpec.b_lpec = full(lpec_functions.b_res_fun(x_lin,p0,M_val));
    n_biactive = dims.n_comp;
else
    A_lpec = full(lpec_functions.A_lpec_fun(M_val));
    % A_lpec = lpec_functions.A_lpec;
    % TODO: check tolerance for active set min(p0_relaxed(end),tol)
    active_set_estimate_k = find_active_sets(x_lin, dims, tol_active_set);
    n_biactive = sum(active_set_estimate_k.I_00);
        lpec.lb(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
    lpec.lb(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
    if settings.reduced_lpec_via_fixed_integers
        lpec.lb_binary(active_set_estimate_k.I_plus_0) = 1;
        lpec.ub_binary(active_set_estimate_k.I_plus_0) = 1;
        lpec.lb_binary(active_set_estimate_k.I_0_plus) = 0;
        lpec.ub_binary(active_set_estimate_k.I_0_plus) = 0;
        b_lpec = full(lpec_functions.b_res_fun(x_lin,p0,M_val));
    else
        error('not implemented fully- please set settings.reduced_lpec_via_fixed_integers = true')
        dims.n_auxiliary = n_biactive;
        b_lpec = [x_lin(ind_x1(active_set_estimate_k.I_00));x_lin(ind_x2(active_set_estimate_k.I_00))-M];
        lb_binary_k = 0*ones(n_biactive,1);
        ub_binary_k = 1*ones(n_biactive,1);
        sense_B = repmat('>',1,2*n_biactive);
        sense = [repmat('=',1,dims.n_eq), repmat('>',1,dims.n_ineq), sense_B];
        vtype = [repmat('C',1,n_primal), repmat('B',1,n_biactive)];
        vtype_num = [repmat(0,1,n_primal), repmat(1,1,n_biactive)];
        if n_biactive>0
            A_lpec_k = [A_lpec([find(active_set_estimate_k.I_00);active_set_estimate_k.I_00+n_comp],1:n_primal),A_lpec([1:n_biactive;n_comp+1:n_comp+n_biactive],1:n_biactive)];
        else
            A_lpec_k  = [];
        end
    end
    lpec.A_lpec = A_lpec;
    lpec.b_lpec = b_lpec;
end


% generic initalization, these three are always updated in the loop:
lpec.d_lpec = zeros(dims.n_primal,1); % initial guess for cont. variables
lpec.y_lpec = zeros(dims.n_comp,1); % inital guess for bin. variablels.
lpec.rho_TR = rho_TR;

% lpec.f = nabla_f; % cost gradient
% lpec.A_eq = A_eq_k;
% lpec.b_eq = b_eq_k;
% lpec.A_ineq = A_ineq_k;
% lpec.b_ineq = b_ineq_k;
% lpec.A_lpec = A_lpec_k; % constraint matrices for the binary constraint to model complementarities
% lpec.b_lpec = b_lpec_k;
end