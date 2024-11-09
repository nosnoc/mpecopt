function [results_lpec,stats_lpec] = lpec_fallback_strategy(lpec,settings,results_lpec,stats_lpec)

% TODO: Add verbose for residuals
settings_lpec = settings.settings_lpec;
ii = 0;
solver_message_in =  stats_lpec.solver_message;
solver_message_out = solver_message_in;
if settings.verbose_solver
    print_iter_line();
end
% initial guess
if ~any(isnan(results_lpec.y_lpec))
    lpec.y_lpec = results_lpec.y_lpec;
end
if ~any(isnan(results_lpec.d_lpec))
    lpec.d_lpec = results_lpec.d_lpec;
end
n_lpecs_solved = 0;

stopping_criterion_satisfied = false;
if strcmp(solver_message_in,'INFEASIBLE')
    while ii < settings.N_attemps_infeasible && ~stopping_criterion_satisfied 
        
        if ii ~= settings.N_attemps_infeasible-1
            lpec.rho_TR = settings.fallback_TR_increasing_factor*lpec.rho_TR;
        else
            lpec.rho_TR = max(settings.fallback_TR_increasing_factor*lpec.rho_TR,5*max(abs(lpec.x_lin)));
            lpec.b_lpec(abs(lpec.b_lpec) > settings.BigM/2) = lpec.b_lpec(abs(lpec.b_lpec) > settings.BigM/2)*1e4;
        end
        [results_lpec,stats_lpec] = lpec_solver(lpec,settings_lpec);
        n_lpecs_solved = n_lpecs_solved+1;
        ii = ii+1;
        if settings.verbose_solver
            print_iter_stats('R',ii,results_lpec.f_opt,0,0,'LPEC-R-IN',stats_lpec.nodecount,stats_lpec.solver_message,lpec.rho_TR,norm(results_lpec.d_lpec),stats_lpec.lpec_solution_exists,0)
        end
        % warm start
        lpec.d_lpec = results_lpec.d_lpec;
        lpec.y_lpec = results_lpec.y_lpec;
        solver_message_out = stats_lpec.solver_message;
        % if strcmp(solver_message_out,'OPTIMAL') || (strcmp(solver_message_out,'NODE_LIMIT') && isnan(results_lpec.f_opt))
        if ~isnan(results_lpec.f_opt)
            stopping_criterion_satisfied  = true;
        end
    end
end

if strcmp(solver_message_in,'NODE_LIMIT')
    F_iters = [results_lpec.f_opt];
    D_iters = [results_lpec.d_lpec];
    Y_iters = [results_lpec.y_lpec];
    TR_iters = [lpec.rho_TR];

    while ii <= settings.N_attemps_node_limit && ~stopping_criterion_satisfied 
        lpec.rho_TR = settings.fallback_TR_reducing_factor*lpec.rho_TR;
        [results_lpec,stats_lpec] = lpec_solver(lpec,settings_lpec);
        n_lpecs_solved = n_lpecs_solved+1;
        if settings.verbose_solver
            try
             print_iter_stats('R',ii,results_lpec.f_opt,0,0,'LPEC-R-NL',stats_lpec.nodecount,stats_lpec.solver_message,lpec.rho_TR,norm(results_lpec.d_lpec),stats_lpec.lpec_solution_exists)
            catch
                fprintf("verbose failed\n");
            end
        end
        ii = ii+1;
        % warm start
        lpec.d_lpec = results_lpec.d_lpec;
        lpec.y_lpec = results_lpec.y_lpec;
        solver_message_out = stats_lpec.solver_message;
        F_iters = [F_iters,results_lpec.f_opt];
        D_iters = [D_iters,results_lpec.d_lpec];
        Y_iters = [Y_iters,results_lpec.y_lpec];
        TR_iters = [TR_iters, lpec.rho_TR ];
        if strcmp(solver_message_out,'OPTIMAL') 
            % || (strcmp(solver_message_out,'NODE_LIMIT') && isnan(results_lpec.f_opt))
            stopping_criterion_satisfied  = true;
        end
    end
    [val,ind] = min(F_iters);
    results_lpec.f_lpec = F_iters(:,ii);
    results_lpec.d_lpec = D_iters(:,ii);
    results_lpec.y_lpec = Y_iters(:,ii);
    lpec.rho_TR = TR_iters(:,ii);
    % a

end

if settings.verbose_solver
    print_iter_line();
end
results_lpec.n_lpecs_solved = n_lpecs_solved;
results_lpec.rho_TR = lpec.rho_TR;
end
