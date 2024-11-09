if settings.use_lpec_fallback_strategy_phase_ii && ismember(settings.settings_lpec.lpec_solver,{'Gurobi','Highs','Matlab'})
    if strcmp(stats_lpec.solver_message,'INFEASIBLE') || (strcmp(stats_lpec.solver_message,'NODE_LIMIT') && isnan(results_lpec.f_opt))
        if settings.debug_mode_on
            keyboard;
        end
        [results_lpec,stats_lpec] = lpec_fallback_strategy(lpec,settings,results_lpec,stats_lpec);
        n_lpec_total = n_lpec_total + results_lpec.n_lpecs_solved; n_lpec_k = n_lpec_k + results_lpec.n_lpecs_solved;
        lpec_solution_exists = stats_lpec.lpec_solution_exists;
        d_lpec_k_l = results_lpec.d_lpec;  y_lpec_k_l = results_lpec.y_lpec; f_lin_opt_k_l = results_lpec.f_opt;
        x_trail_lpec = x_k + d_lpec_k_l;
        h_comp_lpec_k_l = full(h_comp_fun(x_trail_lpec,p0));
        h_std_lpec_k_l = full(h_std_fun(x_trail_lpec,p0));
        if settings.take_TR_from_fallback_strategy
            rho_TR_k_l = results_lpec.rho_TR;
        end
    end
end