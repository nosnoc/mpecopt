if settings.fast_B_stationarity_check
    % && l_k == 1
    % If this is active, solve a LPEC with a small TR, if B-stat stop, else use it as initla guess
    [lpec_is_feasible,min_rho_TR,lpec_infeasiblity] = check_feasiblity_of_lpec(lpec,settings); % check if d = 0 is feasible;
    if lpec_is_feasible
        lpec.d_lpec = d_lpec_k_l;  lpec.y_lpec = y_lpec_k_l; % inital guess for bin. variablels.
        lpec.rho_TR = 10*settings.tol_B_stationarity; % update trust region - small TR to isolate d = 0;
        [results_lpec,stats_lpec] = lpec_solver(lpec,settings.settings_lpec);
        d_lpec_k_l = results_lpec.d_lpec;  y_lpec_k_l = results_lpec.y_lpec; f_lin_opt_k_l = results_lpec.f_opt;
        h_total_k = full(h_total_fun(x_k,p0));
        [success,stopping_criterion_fullfiled,solver_message] = check_if_B_stationary(d_lpec_k_l,f_lin_opt_k_l,nabla_f_k,h_total_k,success,stopping_criterion_fullfiled,solver_message,settings);
        if success
            b_stationarity = true;
        end
    end
end