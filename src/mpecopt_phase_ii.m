function [solution,stats] = mpecopt_phase_ii(mpec_casadi,lpec_casadi,dims,settings,solver_initalization,stats,phase_ii)
% This function implements the Phase II algorithm of MPECppt.
k = 1;
n_cycles = 0;
resolve_nlp = true;
x_k = solver_initalization.x0;
unfold_struct(mpec_casadi,'caller');
unfold_struct(solver_initalization,'caller');
% unfold_struct(mpec_casadi,'caller');
if phase_ii
    rho_TR_k_l = settings.rho_TR_phase_ii_init;
else
    rho_TR_k_l = settings.rho_TR_phase_i_init;
    rho_TR_k_l = max(rho_TR_k_l,2*max(abs(x_k)));
    if full(h_comp_fun(x_k,p0)) > settings.tol_feasibility && ~settings.feasibility_project_to_bounds
        % Make sure that Phase I applied to feasbility MPECs has TR a large enough trust region to satisfiy inital comp. constraints;
        rho_TR_init_lb = 2*max(abs(min((x_k(dims.ind_x1)),(x_k(dims.ind_x2)))));
        % rho_TR_init_ub = 2*max(max(abs(x_k(dims.ind_x1)),abs(x_k(dims.ind_x2))));
        rho_TR_k_l = max(rho_TR_k_l,rho_TR_init_lb);
    end
end
% Check if problem infeasible, then add proper TR that allows feasiblity;
if phase_ii
t_phase_ii_start = tic;
end
% --------------------------- major (outer) itrations -------------------------------------------
while (k <= settings.max_iter) && ~stats.stopping_criterion_fullfiled
    l_k = 1;  % minor/inner iter counter in k-th major/outer iteration
    n_nlp_k = 0; n_lpec_k = 0; % lpecs/nlps solved in iteration k
    accept_trail_step = false;
    if settings.reset_TR_radius && k > 1 
        if phase_ii
            rho_TR_k_l = settings.rho_TR_phase_ii_init;
        else
            rho_TR_k_l = settings.rho_TR_phase_i_init;
        end
    else
        rho_TR_k_l = min(rho_TR_k_l, settings.rho_TR_max); % initalize TR for inner loop
    end

    f_k  = full(mpec_casadi.f_fun(x_k,p0));
    h_comp_k_l = full(mpec_casadi.h_comp_fun(x_k,p0));
    h_std_k_l = full(mpec_casadi.h_std_fun(x_k,p0));
    f_opt_k_l = full(mpec_casadi.f_fun(x_k,p0));
    nabla_f_k  = full(mpec_casadi.nabla_f_fun(x_k,p0));
    if norm(nabla_f_k) >= 1e2 && settings.rescale_large_objective_gradients
        nabla_f_k = nabla_f_k./(norm(nabla_f_k));
        lpec.f = nabla_f_k;
    end
    t_lpec_preparation_iter = tic;
    lpec = create_lpec_subproblem(x_k,p0,rho_TR_k_l,lpec_casadi,dims,settings,settings.tol_active);
    stats.iter.cpu_time_lpec_preparation_iter = [stats.iter.cpu_time_lpec_preparation_iter;toc(t_lpec_preparation_iter)];

    % --------------------------- Inner (minor) itrations ------------------------------------------
    while ~accept_trail_step && l_k <= settings.max_inner_iter && ~stats.stopping_criterion_fullfiled
        % Here one could do a fast_B_stationarity_check
        y_lpec_k_previous = y_lpec_k_l; % to keep track of active set chnages
        %lpec.d_lpec = d_lpec_k_l; % Initial guess and TR for the LPEC
        lpec.y_lpec = y_lpec_k_l; % inital guess for bin. variablels.
        lpec.rho_TR = rho_TR_k_l; % update trust region
        stats.iter.rho_TR_iter = [stats.iter.rho_TR_iter, rho_TR_k_l]; % store TR radius
        % Solve LPEC
        [results_lpec,stats_lpec] = lpec_solver(lpec,settings.settings_lpec);
        if ~phase_ii
            stats.iter.cpu_time_lpec_phase_i_iter = [stats.iter.cpu_time_lpec_phase_i_iter, stats_lpec.cpu_time]; % stats
        else
            stats.iter.cpu_time_lpec_phase_ii_iter = [stats.iter.cpu_time_lpec_phase_ii_iter, stats_lpec.cpu_time]; % stats
        end
        n_lpec_k = n_lpec_k + 1; stats.n_lpec_total = stats.n_lpec_total + 1;
        % extract LPEC results
        lpec_solution_exists = stats_lpec.lpec_solution_exists;
        d_lpec_k_l = results_lpec.d_lpec;
        y_lpec_k_l = results_lpec.y_lpec;
        f_lin_opt_k_l = results_lpec.f_opt;
        stats.f_lpec = results_lpec.f_opt; 
        x_trail_lpec = x_k + d_lpec_k_l;
        % Infeasiblity check
        h_comp_lpec_k_l = full(mpec_casadi.h_comp_fun(x_trail_lpec,p0));
        h_std_lpec_k_l = full(mpec_casadi.h_std_fun(x_trail_lpec,p0));
        if settings.verbose_solver
            print_iter_stats(k,l_k,f_lin_opt_k_l,h_std_lpec_k_l,h_comp_lpec_k_l,'LPEC',stats_lpec.nodecount,stats_lpec.solver_message,lpec.rho_TR,norm(d_lpec_k_l),stats_lpec.cpu_time,' ')
        end
        % One could do a fallbakc strategy here if the LPEC fails, but should not happen as these LPECs are always feasible, code can be found in lpec_fallback_strategy_phase_ii.m
        if lpec_solution_exists
            if settings.plot_lpec_iterate
                plot_lpec(nabla_f_k, x_k, d_lpec_k_l, rho_TR_k_l)
            end
            % --------------------------- Check if B-stationary point found --------------------------
            h_total_k = full(mpec_casadi.h_total_fun(x_k,p0));
            if (h_total_k <= settings.tol) && ((abs(f_lin_opt_k_l) <= settings.tol_B_stationarity || norm(nabla_f_k) <= settings.tol_B_stationarity))  % if objective zero (either if cost gradient zero, or solution leads to it) = then set step to zero => B stationarity
                if settings.reset_lpec_objective
                    d_lpec_k_l = d_lpec_k_l*0; % if the current point is feasible, and the objective is zero, then d = 0 is also a solution of the lpec (occurs if a solution is not on the verties of the lp)
                    f_lin_opt_k_l = 0;
                end
            end
            if norm(d_lpec_k_l) <= settings.tol_B_stationarity
            % if abs(f_lin_opt_k_l) <= settings.tol_B_stationarity 
                stats.stopping_criterion_fullfiled = true;     % B-stationary point found, optimal solution found!
                stats.solver_message = 'B-stationary point found sucessfully.';
                stats.success = true;
                stats.b_stationarity = true;
                resolve_nlp = false;
                if settings.count_first_lpec_into_phase_i && k == 1 && l_k == 1 && phase_ii
                    stats.solved_in_phase_i = true;
                end
            end

            stats.iter.X_lpec = [stats.iter.X_lpec, x_trail_lpec];
            stats.iter.d_norm_lpec = [stats.iter.d_norm_lpec, norm(d_lpec_k_l)];
            stats.iter.f_lpec = [stats.iter.f_lpec, f_lin_opt_k_l]; % store some stats
            % --------------------------- set up piece NLP with new active set-------------------------
            if settings.compute_bnlp_step && ~stats.stopping_criterion_fullfiled
                lbx_bnlp_k = lbx; ubx_bnlp_k = ubx;  % reset bounds of bnlp.
                lbg_tnlp_k = lbg; ubg_tnlp_k = ubg;
                % find_active_sets_tnlp
                active_set_estimate_k = find_active_sets_piece_nlp(x_trail_lpec,nabla_f_k,y_lpec_k_l,dims,settings,settings.tol_active);
                ubx_bnlp_k(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
                ubx_bnlp_k(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
                ubx_bnlp_k(dims.ind_x1(active_set_estimate_k.I_00)) = 0;
                ubx_bnlp_k(dims.ind_x2(active_set_estimate_k.I_00)) = 0;
                active_set_changes_k_l = round(sum(abs(y_lpec_k_l-y_lpec_k_previous)));
                stats.iter.active_set_changes  = [stats.iter.active_set_changes, active_set_changes_k_l];

                if active_set_changes_k_l > 0 || (l_k == 1)
                    resolve_nlp  = true;
                else
                    resolve_nlp  = false;
                end
                % --------------------------- solve piece NLP -------------------------
                if resolve_nlp
                    t_nlp_start = tic;
                    results_nlp = solver('x0',x_k,'p', p0, 'lbx', lbx_bnlp_k, 'ubx', ubx_bnlp_k,'lbg', lbg_tnlp_k, 'ubg', ubg_tnlp_k);
                    cpu_time_nlp_k_l = toc(t_nlp_start);
                    x_trail_nlp = full(results_nlp.x);
                    lambda_x_trail_nlp = full(results_nlp.lam_x);
                    stats_nlp = solver.stats();
                    nlp_iters_k_l = stats_nlp.iter_count;
                    h_comp_k_l = full(mpec_casadi.h_comp_fun(x_trail_nlp,p0));
                    h_std_k_l = full(mpec_casadi.h_std_fun(x_trail_nlp,p0));
                    f_opt_k_l = full(mpec_casadi.f_fun(x_trail_nlp,p0));
                    f_k_trail = full(mpec_casadi.f_fun(x_trail_nlp, p0));
                    n_nlp_k = n_nlp_k+1;
                    d_nlp_k_l = x_trail_nlp-x_k;
                    stats.n_nlp_total = stats.n_nlp_total+1;
                    if ~phase_ii
                        stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter, cpu_time_nlp_k_l];
                    else
                        stats.iter.cpu_time_nlp_phase_ii_iter = [stats.iter.cpu_time_nlp_phase_ii_iter, cpu_time_nlp_k_l];
                    end
                    stats.iter.X_all = [stats.iter.X_all, x_trail_nlp];
                    nlp_step_sucessful = true;

                    if isequal(stats_nlp.return_status,'Infeasible_Problem_Detected') && settings.debug_mode_on 
                        % this may hapen if lpec makes big jump and xk not fesaible for new bnlp
                        keyboard;                   
                    end

                    if ~ismember(stats_nlp.return_status, {'Solve_Succeeded', 'Search_Direction_Becomes_Too_Small', 'Solved_To_Acceptable_Level'})
                    % if full(mpec_casadi.h_total_fun(x_trail_nlp,p0)) <= settings.tol_feasibility
                        % Remark: This may happen if the full LPEC is allowed to jump over inactive corners, otherwise, with the reduced LPEC it always starts from a feasible point.
                        % stopping_criterion_fullfiled = true;
                        accept_trail_step = false;
                        nlp_step_sucessful = false;
                        solver_message = ['NLP solver faild to solve BNLP in Phase II: ' stats_nlp.return_status '.'];
                        if settings.debug_mode_on
                            keyboard;
                        end
                    else
                        % -------------------- Check is the BNLP solution S-stationary ------------------------------
                        if settings.stop_if_S_stationary && ~stats.stopping_criterion_fullfiled
                            active_set_estimate_k = find_active_sets(x_trail_nlp, dims, settings.tol_active); % check is it S-stationary;
                            if sum(active_set_estimate_k.I_00) == 0
                                stats.stopping_criterion_fullfiled = true;
                                stats.multiplier_based_stationarity = 'S';
                                stats.solver_message = 'S-stationary point found, BNLP biactive set is empty.';
                                stats.b_stationarity = true;
                                stats.success = true;
                                stats.f_lpec = 0;
                                x_k = x_trail_nlp;
                            else
                                if settings.resolve_tnlp_for_stationarity && ~settings.PieceNLPStartegy == "TNLP"
                                    % get rid of bilinear constraints
                                    ubg_relaxed(ind_comp) = +inf;
                                    if strcmp(settings.relax_and_project_homotopy_parameter_steering,'Ell_1') || strcmp(settings.relax_and_project_homotopy_parameter_steering,'Ell_inf')
                                        ubx_relaxed(end) = +inf;
                                    end
                                    ubx_relaxed(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
                                    ubx_relaxed(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
                                    ubx_relaxed(dims.ind_x1(active_set_estimate_k.I_00)) = 0;
                                    ubx_relaxed(dims.ind_x2(active_set_estimate_k.I_00)) = 0;
                                    % solve TNLP
                                    results_nlp = solver_relaxed('x0',x_k_init,'p',p0_relaxed,'lbx',lbx_relaxed,'ubx',ubx_relaxed,'lbg',lbg_relaxed,'ubg',ubg_relaxed);
                                    x_k = full(results_nlp.x);
                                    x_k = x_k(1:dims.n_primal);
                                    lambda_x_k = full(results_nlp.lam_x);
                                    lambda_x_k  = lambda_x_k(1:dims.n_primal);
                                    % lambda_g_k = full(results_nlp.lam_g);
                                    [stats.multiplier_based_stationarity, solver_message] = determine_multipliers_based_stationary_point(x_k,lambda_x_k,dims,settings);
                                    if strcmp(stats.multiplier_based_stationarity,'S')
                                        stats.stopping_criterion_fullfiled = true;
                                        stats.solver_message = 'S-stationary point found, BNLP biactive set is nonempty.';
                                        stats.success = true;
                                        stats.b_stationarity = true;
                                        stats.f_lpec = 0;
                                    end
                                end
                            end
                        end
                    end
                else
                    d_nlp_k_l = 0*d_nlp_k_l;
                end
            else
                if ~settings.compute_bnlp_step
                    x_trail_nlp = x_trail_lpec;  % take only lpec steps (slower convergence) % only linear steps (needs apropiate globalization strategy)
                    f_k_trail = full(mpec_casadi.f_fun(x_trail_nlp, p0));
                    stats.iter.X_all = [stats.iter.X_all, x_trail_nlp];
                end
            end
            % --------------------------- Globalization: check step acceptence -------------------------
            t_globalization_iter = tic;
            if ~stats.stopping_criterion_fullfiled 
                % Current iterate
                f_k = full(mpec_casadi.f_fun(x_k, p0));
                h_std_k = full(mpec_casadi.h_std_fun(x_k, p0));
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k, p0));
                h_total_k = max(h_std_k,h_comp_k);
                % Trail point
                f_k_trail = full(mpec_casadi.f_fun(x_trail_nlp, p0));
                delta_f_k_l = f_k-f_k_trail;
                h_std_k_trail = full(mpec_casadi.h_std_fun(x_trail_nlp, p0));
                h_comp_k_trail = full(mpec_casadi.h_comp_fun(x_trail_nlp, p0));
                h_total_k_trail = max(h_std_k_trail,h_comp_k_trail);
                % f_k_linear_trail = full(mpec_casadi.f_fun(x_k+d_lpec_k_l, p0));

                % Check if sufficient objective decerase w.r.t to current iterate
                feasible_enough = h_total_k_trail < settings.tol_feasibility;
                if feasible_enough
                    if f_k_trail < f_k
                        accept_trail_step = true;
                    else
                        accept_trail_step = false;
                        if settings.debug_mode_on
                            % debug rejected step;
                            keyboard;
                        end
                    end
                else
                    % in case Phase I gave a point that is not feasible, but phase II fixed it - normally this should never happen.
                    if nlp_step_sucessful
                        accept_trail_step = true;
                    end
                    if settings.debug_mode_on
                        keyboard;
                    end
                end
                
                % Update TR:
                if accept_trail_step 
                    rho_TR_k_l = settings.TR_increasing_factor*rho_TR_k_l;
                else
                    rho_TR_k_l = settings.TR_reducing_factor*rho_TR_k_l;
                end
                % rho_TR_k_l = max(settings.rho_TR_min,rho_TR_k_l);
                if rho_TR_k_l < settings.rho_TR_min && l_k > 1
                    n_cycles = n_cycles+1;
                    break; % avoid cyclin with rho_TR_min
                end
                
                stats.iter.cpu_time_globalization_iter = [stats.iter.cpu_time_globalization_iter; toc(t_globalization_iter)];
                % ------------------- Debug mode - catch some unexpected behaviour ---------------------
                if ~feasible_enough && k>1 && settings.debug_mode_on
                    keyboard; % Stop if Phase II has infeasible problems;
                end
                if feasible_enough && f_k_trail > f_k && settings.debug_mode_on
                    % Stop here if the objective incerased when started from a feasible point, see also % results_lpec.f_opt;
                    % This can happen if LPEC predictes reduction, but actual redicution does not happen
                    keyboard;
                end
                % reject step, reduce TR and compute new step (check the filte only if nlp resolved)
                if settings.accept_last_inner_iter && l_k == settings.max_inner_iter && ~accept_trail_step
                    accept_trail_step = true;
                    if settings.debug_mode_on
                        keyboard;
                    end
                end
            end
            if accept_trail_step
                x_k = x_trail_nlp; % completion of one outer iter
                if settings.compute_bnlp_step
                    lambda_x_k = full(results_nlp.lam_x);
                end
            end
            % ----------------------------  verbose current inner NLP iter ------------------------------
            if settings.verbose_solver && resolve_nlp && settings.compute_bnlp_step
                % print_iter_stats(k,l_k,f_opt_k_l,h_std_k_l,h_comp_k_l,'BNLP',nlp_iters_k_l,stats_nlp.return_status,rho_TR_k_l,norm(d_nlp_k_l),stats.multiplier_based_stationarity,accept_trail_step) % Verbose inner iterations.
                print_iter_stats(k,l_k,f_opt_k_l,h_std_k_l,h_comp_k_l,'BNLP',nlp_iters_k_l,stats_nlp.return_status,rho_TR_k_l,delta_f_k_l,cpu_time_nlp_k_l,accept_trail_step) % Verbose inner iterations.
            end
        else
            % LPEC failed -- This can happen only if a CQ does not hold or an infeasible point was passed to Phase II
            stats.solver_message = 'LPEC inconsistent - solution failed.';
            stats.problem_infeasible = true;
            stats.stopping_criterion_fullfiled = true;

            if settings.debug_mode_on
                keyboard;
            end
        end
        l_k = l_k+1; % inner iterations counter
    end
    % compute new objective and residuals - with accepted point
    f_k = full(mpec_casadi.f_fun(x_k,p0));
    h_std_k = full(mpec_casadi.h_std_fun(x_k,p0));
    h_comp_k = full(mpec_casadi.h_comp_fun(x_k,p0));
    % store stats
    stats.total_inner_iterations = [stats.total_inner_iterations, l_k];
    stats.iter.X_outer  = [stats.iter.X_outer, x_k];
    stats.iter.f_iter = [stats.iter.f_iter, f_k];
    stats.iter.h_std_iter = [stats.iter.h_std_iter, h_std_k];
    stats.iter.h_comp_iter = [stats.iter.h_comp_iter, h_comp_k];
    stats.iter.n_nlp_iter = [stats.iter.n_nlp_iter, n_nlp_k];
    stats.iter.n_lpec_iter = [stats.iter.n_lpec_iter, n_lpec_k];
    k = k+1;
    if settings.verbose_solver
        print_iter_line();
    end
    if n_cycles >= 3
        stats.solver_message = 'Major loop was cycling due to bad problem scaling or too low tolerances.';
        break;
    end
end
stats.rho_TR_final = rho_TR_k_l;
if phase_ii
stats.cpu_time_phase_ii = toc(t_phase_ii_start);
end

if k >= settings.max_iter
    stats.max_iterations_reached = true;
    if settings.debug_mode_on
        keyboard;
    end
end
%  ------------- max iteration but early terminaton tolorance achieved?---------------------
if ((stats.max_iterations_reached && stats.success == 0) || n_cycles == 3)&& settings.allow_early_termination
    if (h_total_k <= settings.tol_B_stationarity_early_term) && ((abs(f_lin_opt_k_l) <= settings.tol_B_stationarity_early_term|| norm(nabla_f_k) <= settings.tol_B_stationarity_early_term))  % if objective zero (either if cost gradient zero, or solution leads to it) = then set step to zero => B stationarity
        % B-stationary point found, optimal solution found!
        if n_cycles == 3
            stats.solver_message = 'Major loop was cycling due to bad problem scaling or too low tolerances. B-stationary point found at lower tolerance.';
        else
            stats.solver_message = 'Maximum number of iterations reached, B-stationary point found at lower tolerance.';
        end
        stats.success = true;
        stats.b_stationarity = true;
        stats.f_lpec = f_lin_opt_k_l;
    end
    if settings.debug_mode_on
        keyboard;
    end
end

% --------------- compute multiplier-based stationary points --------------
multiplier_based_stationarity_debug = stats.multiplier_based_stationarity;
if (stats.success || k==settings.max_iter) && settings.compute_tnlp_stationary_point && phase_ii
    % if ~strcmp(settings.piece_nlp_strategy,'TNLP')
        % resolve TNLP for correct multipliers
        lbx_bnlp_k = lbx; ubx_bnlp_k = ubx;  % reset bounds of bnlp.
        lbg_tnlp_k = lbg; ubg_tnlp_k = ubg;
        initial_strategy = settings.piece_nlp_strategy;
        settings.piece_nlp_strategy = 'TNLP'; % used to get tnlp active sets
        nabla_f_k = full(mpec_casadi.nabla_f_fun(x_k,p0));
        active_set_estimate_k = find_active_sets_piece_nlp(x_k,nabla_f_k,y_lpec_k_l,dims,settings,settings.tol_active);
        ubx_bnlp_k(dims.ind_x1(active_set_estimate_k.I_0_plus)) = 0;
        ubx_bnlp_k(dims.ind_x2(active_set_estimate_k.I_plus_0)) = 0;
        ubx_bnlp_k(dims.ind_x1(active_set_estimate_k.I_00)) = 0;
        ubx_bnlp_k(dims.ind_x2(active_set_estimate_k.I_00)) = 0;
        settings.piece_nlp_strategy = initial_strategy;
    % end
        t_nlp_start = tic;
        results_nlp = solver('x0',x_k,'p', p0, 'lbx', lbx_bnlp_k, 'ubx', ubx_bnlp_k,'lbg', lbg_tnlp_k, 'ubg', ubg_tnlp_k);
        cpu_time_nlp_k_l = toc(t_nlp_start);
        x_k_multi = full(results_nlp.x);
        lambda_x_k = full(results_nlp.lam_x);
        [stats.multiplier_based_stationarity, ~] = determine_multipliers_based_stationary_point(x_k_multi,lambda_x_k,dims,settings);
end

% Debug falure of stationary point computation
if ~strcmp(multiplier_based_stationarity_debug,'X') && settings.debug_mode_on &&  ~isequal(multiplier_based_stationarity_debug,stats.multiplier_based_stationarity)
    keyboard;
end

% eval stats at solution
% x_k = x_k(1:dims.n_primal_non_lifted);

solution.x = x_k(1:dims.n_primal_non_lifted);
solution.x_lifted = x_k;
solution.f = full(mpec_casadi.f_fun(x_k,p0));
solution.x1_opt = x_k(dims.ind_x1);
solution.x2_opt = x_k(dims.ind_x2);
stats.h_std = full(mpec_casadi.h_std_fun(x_k,p0));
stats.comp_res = full(mpec_casadi.h_comp_fun(x_k,p0));
stats.h_total = full(mpec_casadi.h_total_fun(x_k,p0));
stats.total_outer_iterations = k;

stats.n_biactive = sum(x_k(dims.ind_x1)+x_k(dims.ind_x2) < 2*settings.tol_active );
try
    stats.lambda_g_opt = full(results_nlp.lam_g);
    stats.lambda_x_opt = full(results_nlp.lam_x);
    stats.n_active_ineq = sum(abs(lambda_g_opt(ind_g_ineq))>settings.tol);
    stats.n_active_box = sum(abs((lambda_x_opt))>settings.tol & lbx_bnlp_k~=ubx_bnlp_k);
    stats.n_box = sum(lbx_bnlp_k~=ubx_bnlp_k);
    stats.n_box_simple = sum(lbx_bnlp_k==ubx_bnlp_k);
catch
    stats.n_active_ineq = nan;
    stats.n_active_box = nan;
    stats.n_box = sum(lbx~=ubx);
    stats.n_box_simple = sum(lbx==ubx);
end

end

