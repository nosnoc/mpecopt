classdef Mpecopt < handle
    properties
        settings
        mpec
        mpec_casadi
        lpec_casadi
        dims
        solver_initialization
        stats
    end

    methods
        function obj = Mpecopt(mpec, solver_initialization, settings)
            obj.settings = settings;
            obj.mpec = mpec;
            [mpec_casadi, dims, solver_initialization, stats] =  create_mpec_functions(mpec,solver_initialization,settings);
            obj.mpec_casadi = mpec_casadi;
            obj.solver_initialization = solver_initialization;
            obj.dims = dims;
            obj.stats = stats;
        end

        function [solution,stats] = solve(obj, solver_initialization)
            arguments
                obj
                solver_initialization = []
            end
            %% Get data
            if isempty(solver_initialization)
                solver_initialization = obj.solver_initialization;
            else
                warning("solver initialization passed, but not processed")
                % TODO(@anton) separate initialization update
                obj.solver_initialization = solver_initialization;
            end
            mpec = obj.mpec;
            settings = obj.settings;
            mpec = obj.mpec;
            mpec_casadi = obj.mpec_casadi;
            lpec_casadi = obj.lpec.casadi;
            dims = obj.dims;
            stats = obj.stats;
            
            %% Ininitialization
            k = 1;
            x_k = solver_initialization.x0;
            y_lpec_k_l = double(x_k(dims.ind_x1)> x_k(dims.ind_x2)); % binary variables
            stats.stopping_criterion_fullfiled = false;
            d_nlp_k_l = x_k;
            d_lpec_k_l = x_k;
            x_trail_nlp = x_k;
            feasible_bnlp_found = false; % did phase i find a feasible bnlp?

            rho_TR_k_l = settings.rho_TR_phase_i_init;
            % stats.iter.X_outer = x_k; % save all accepted steps;
            % stats.iter.X_all = x_k;  % save all iterates of accepted and rejceted steps;
            stats.iter.X_outer = []; % save all accepted steps;
            stats.iter.X_all = [];  % save all iterates of accepted and rejceted steps;
            stats.iter.X_lpec = []; % save all lpec trail points
            stats.iter.d_norm_lpec = []; % trail point norms
            stats.iter.f_lpec = []; % objectives of the lpecs
            stats.total_inner_iterations = [];
            stats.iter.rho_TR_iter = [];
            stats.iter.f_iter = []; % objective of accepted steps
            stats.iter.h_std_iter = []; % standard infeasilibty
            stats.iter.h_comp_iter = []; % comp infeasilibty
            stats.n_nlp_total = 0;
            stats.n_lpec_total = 0;
            stats.iter.n_nlp_iter = [];
            stats.iter.n_lpec_iter = [];
            stats.iter.active_set_changes = [];
            stats.solver_message  = 'Maximum number of iterations reached.';
            stats.success = false;
            stats.success_phase_i = false;
            stats.problem_infeasible = false;
            stats.solved_in_phase_i = false;
            stats.multiplier_based_stationarity = 'X'; % not yet determined.
            stats.b_stationarity = false;
            stats.f_lpec = 1;
            stats.rho_TR_final = 1e-3;
            stats.max_iterations_reached = false;
            resolve_nlp = true; % after every LPEC solve in the main loop, solve new NLP;
            rho_TR_k_l = 1e-3;

            if settings.plot_mpec_multipliers
                figure;     % TODO: add figure handle that is passed to the multipl function.
            end
            cpu_time_prepare_mpec = toc(t_prepare_mpec);

            % ------------------ Prepare homotopy solver for Phase I -----------------------------
            if strcmp(settings.initialization_strategy,"RelaxAndProject")
                t_generate_nlp_solvers = tic;
                % TODO : input just mpec_casadi and solver_initialization
                [solver_relaxed,x_k_init,p0_relaxed,lbx_relaxed,ubx_relaxed,lbg_relaxed,ubg_relaxed] = create_phase_i_nlp_solver(mpec_casadi.f,mpec_casadi.g,mpec_casadi.x,mpec_casadi.x1,mpec_casadi.x2,mpec_casadi.p,lbx,ubx,lbg,ubg,p0,x_k,settings,dims);
                stats.cpu_time_generate_nlp_solvers = stats.cpu_time_generate_nlp_solvers+toc(t_generate_nlp_solvers);
            end

            stats.iter.cpu_time_nlp_phase_i_iter = [0];
            stats.iter.cpu_time_nlp_phase_ii_iter = [0];
            stats.iter.cpu_time_lpec_phase_i_iter = [0];
            stats.iter.cpu_time_lpec_phase_ii_iter = [0];
            stats.iter.cpu_time_lpec_preparation_iter = [];
            stats.iter.cpu_time_nlp_iter = [0];
            stats.iter.cpu_time_globalization_iter = [];

            %% Main piece NLP solver (BNLP or TNLP)
            % t_generate_nlp_solvers = tic;
            % nlp = struct('x', x,'f', f,'g', g,'p',p);
            % solver = nlpsol('solver', 'ipopt', nlp, settings.settings_casadi_nlp);
            % stats.cpu_time_generate_nlp_solvers = toc(t_generate_nlp_solvers);
            % mpec_casadi.solver = solver;

            %%  ---------------------- Phase I - compute feasible point ------------------------
            if settings.verbose_solver
                print_phase_i();
                print_iter_header();
            end
            h_comp_ii = full(mpec_casadi.h_comp_fun(x_k(1:dims.n_primal),p0));
            h_std_ii = full(mpec_casadi.h_std_fun(x_k(1:dims.n_primal),p0));
            f_opt_ii = full(mpec_casadi.f_fun(x_k(1:dims.n_primal),p0));
            print_iter_stats('I',0,f_opt_ii,h_std_ii,h_comp_ii,'/',0,'Initial guess',0,0,0,1)

            t_phase_i_start = tic;
            switch settings.initialization_strategy
              case 'RelaxAndProject'
                h_comp_ii = 1;
                ii = 1;
                y_lpec_k_l = double(x_k(dims.ind_x1)> x_k(dims.ind_x2));
                d_lpec_k_l = zeros(dims.n_primal,1);
                while (ii <= settings.max_recovery_iters) && ~stats.problem_infeasible && ~feasible_bnlp_found
                    switch settings.bnlp_projection_strategy
                      case {'LPEC'}
                        if settings.lpec_solve_if_comp_feasible
                            rho_min = 1.01*max(abs(min(x_k(dims.ind_x1),x_k(dims.ind_x2))));
                            if rho_min > settings.rho_TR_phase_i_init || ii == 1
                                solve_lpec = false;
                            else
                                solve_lpec = true;
                            end
                        else
                            solve_lpec = true;
                        end
                        if solve_lpec
                            % prepare and solve lpec; 
                            t_lpec_preparation_iter = tic;
                            lpec = create_lpec_subproblem(x_k,p0,rho_TR_k_l,lpec_casadi,dims,settings,settings.tol_active);
                            stats.iter.cpu_time_lpec_preparation_iter = [stats.iter.cpu_time_lpec_preparation_iter;toc(t_lpec_preparation_iter)];
                            %  Initial guess and TR for the LPEC
                            y_lpec_k_previous = y_lpec_k_l; % to keep track of active set chnages
                            lpec.d_lpec = d_lpec_k_l;
                            lpec.y_lpec = y_lpec_k_l; % inital guess for bin. variablels.

                            rho_TR_k_l_lb = settings.realx_and_project_scale_factor_rho_tr*max(abs(min(x_k(dims.ind_x1),x_k(dims.ind_x2))));
                            rho_TR_k_l_ub = settings.realx_and_project_scale_factor_rho_tr*max(max(abs(x_k(dims.ind_x1)),abs(x_k(dims.ind_x2))));

                            if settings.relax_and_project_tighter_TR
                                rho_TR_k_l = min(rho_TR_k_l_lb, settings.rho_TR_phase_i_init);
                                rho_TR_k_l = max(1e-4,rho_TR_k_l); % to avoid unerasonably small TR at almost fesaible points;
                                                                   % [rho_TR_k_l_lb,  p0_relaxed(end) ]
                            else
                                rho_TR_k_l = max(rho_TR_k_l_lb, settings.rho_TR_phase_i_init);
                            end
                            lpec.rho_TR  =  rho_TR_k_l;
                            stats.iter.rho_TR_iter = [stats.iter.rho_TR_iter, rho_TR_k_l]; % store TR radius
                                                                                           % Solve LPEC
                            [results_lpec,stats_lpec] = lpec_solver(lpec,settings.settings_lpec);
                            stats.iter.cpu_time_lpec_phase_i_iter = [stats.iter.cpu_time_lpec_phase_i_iter, stats_lpec.cpu_time]; % stats
                            stats.n_lpec_total = stats.n_lpec_total+1;
                            % extract results
                            lpec_solution_exists = stats_lpec.lpec_solution_exists;
                            d_lpec_k_l = results_lpec.d_lpec;
                            y_lpec_k_l = results_lpec.y_lpec;
                            f_lin_opt_k_l = results_lpec.f_opt;
                            stats.f_lpec = f_lin_opt_k_l; 
                            x_trail_lpec = x_k + d_lpec_k_l;
                            % Infeasibility check
                            h_comp_lpec_k_l = full(mpec_casadi.h_comp_fun(x_trail_lpec,p0));
                            h_std_lpec_k_l = full(mpec_casadi.h_std_fun(x_trail_lpec,p0));
                            if settings.verbose_solver
                                print_iter_stats(1,ii,f_lin_opt_k_l,h_std_lpec_k_l,h_comp_lpec_k_l,'LPEC',stats_lpec.nodecount,stats_lpec.solver_message,lpec.rho_TR,norm(d_lpec_k_l),stats_lpec.cpu_time,' ')
                            end
                            if settings.use_lpec_fallback_strategy_phase_i  && ii > ceil(settings.lpec_recovery_start*settings.max_recovery_iters)
                                if strcmp(stats_lpec.solver_message,'INFEASIBLE') || (strcmp(stats_lpec.solver_message,'NODE_LIMIT') && isnan(results_lpec.f_opt))
                                    if settings.debug_mode_on
                                        keyboard;
                                    end
                                    [results_lpec,stats_lpec] = lpec_fallback_strategy(lpec,settings,results_lpec,stats_lpec);
                                    stats.n_lpec_total = stats.n_lpec_total + results_lpec.n_lpecs_solved;
                                    lpec_solution_exists = stats_lpec.lpec_solution_exists;
                                    d_lpec_k_l = results_lpec.d_lpec;  y_lpec_k_l = results_lpec.y_lpec; f_lin_opt_k_l = results_lpec.f_opt;
                                    x_trail_lpec = x_k + d_lpec_k_l;
                                end
                            end
                            %
                            if lpec_solution_exists
                                % --------------------------- Check if B-stationary point found --------------------------
                                h_total_k = full(mpec_casadi.h_total_fun(x_k,p0));
                                nabla_f_k = full(mpec_casadi.nabla_f_fun(x_k,p0));
                                if (h_total_k <= settings.tol) && ((abs(f_lin_opt_k_l) <= settings.tol_B_stationarity || norm(nabla_f_k) <= settings.tol_B_stationarity))  % if objective zero (either if cost gradient zero, or solution leads to it) = then set step to zero => B stationarity
                                    if settings.reset_lpec_objective
                                        d_lpec_k_l = d_lpec_k_l*0; % if the current point is feasible, and the objective is zero, then d = 0 is also a solution of the lpec (occurs if a solution is not on the verties of the lp)
                                        f_lin_opt_k_l = 0;
                                    end
                                end
                                if norm(d_lpec_k_l) <= settings.tol_B_stationarity
                                    stats.stopping_criterion_fullfiled = true;     % B-stationary point found, optimal solution found!
                                    stats.solver_message = 'B-stationary point found sucessfully.';
                                    stats.success = true;
                                    stats.b_stationarity = true;
                                    stats.f_lpec = f_lin_opt_k_l;
                                end
                                stats.iter.X_lpec = [stats.iter.X_lpec, x_trail_lpec];
                                stats.iter.d_norm_lpec = [stats.iter.d_norm_lpec, norm(d_lpec_k_l)];
                                stats.iter.f_lpec = [stats.iter.f_lpec, f_lin_opt_k_l]; % store some stats
                                I_0_plus = y_lpec_k_l==0;
                                I_plus_0 = y_lpec_k_l==1;
                                % I_00 = [];
                                active_set_guess_exists = true;
                            else
                                active_set_guess_exists = false;
                            end
                        else
                            active_set_guess_exists = false;
                        end

                      case 'Simple'
                        I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                        I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                        active_set_guess_exists = true;
                      case 'ActiveSetFunction'
                        % TODO: use some fancy function like Lin Fukushima 2006
                        error('not supported at the moment, choos different option');
                        % active_set_identification_fun;
                        % active_set_identification_fun =
                        % tau > 0, eta \in (0,1)
                        tau = 1e-1;eta = 0.5;
                        active_set_identification  = tau*norm(min(x(dims.ind_x1),x(dims.ind_x2)))^eta;
                        active_set_identification_fun = Function('active_set_identification_fun',{x},{active_set_identification});
                        % alpha = I_+0 , beta I_00, gamma I_0+
                        rho_as = full(active_set_identification_fun(x_k));
                        I_0_plus  = x_k(dims.ind_x1) <=  rho_as & x_k(dims.ind_x2) > rho_as ;
                        I_plus_0  = x_k(dims.ind_x1) > rho_as & x_k(dims.ind_x2) <= rho_as ;
                        I_00  = x_k(dims.ind_x1) <= rho_as  & x_k(dims.ind_x2) <= rho_as ;
                        % split I_00
                        I_0_plus((x_k(dims.ind_x1(I_00))>x_k(dims.ind_x2(I_00)))) = 1;
                        I_plus_0((x_k(dims.ind_x1(I_00))<=x_k(dims.ind_x2(I_00)))) = 1;
                        active_set_guess_exists = true;
                    end
                    % solve BNLP/NLP
                    h_total_k = full(mpec_casadi.h_total_fun(x_k,p0));
                    if h_total_k <= settings.tol_feasibility && ii>1
                        feasible_bnlp_found = true;
                    end

                    if active_set_guess_exists && ~feasible_bnlp_found
                        lbx_bnlp_k = lbx;
                        ubx_bnlp_k = ubx;
                        ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                        ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                        t_presolve_nlp_iter = tic;
                        results_nlp = mpec_casadi.solver('x0',x_k,'p',p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',lbg,'ubg',ubg);
                        cpu_time_bnlp_ii = toc(t_presolve_nlp_iter);
                        stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter; cpu_time_bnlp_ii];
                        x_k = full(results_nlp.x);
                        lambda_x_k  = full(results_nlp.lam_x);
                        stats_nlp = mpec_casadi.solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                        h_comp_k = full(mpec_casadi.h_comp_fun(x_k,p0));
                        h_std_k = full(mpec_casadi.h_std_fun(x_k,p0));
                        f_opt_k = full(results_nlp.f);
                        if settings.verbose_solver
                            print_iter_stats(1,ii,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,rho_TR_k_l,norm(x_k_init(1:dims.n_primal)-x_k),cpu_time_bnlp_ii ,1)
                        end
                        if full(mpec_casadi.h_total_fun(x_k,p0)) <= settings.tol_feasibility
                            % if ismember(stats_nlp.return_status, {'Solve_Succeeded', 'Search_Direction_Becomes_Too_Small', 'Solved_To_Acceptable_Level'})
                            feasible_bnlp_found = true;
                        else
                            feasible_bnlp_found = false;
                        end
                    end
                    % If not sucessful do one more relaxation step until max reached
                    if feasible_bnlp_found
                        stats.success_phase_i = true;
                        break;
                    else
                        t_presolve_nlp_iter = tic;
                        results_nlp = solver_relaxed('x0',x_k_init,'p',p0_relaxed,'lbx',lbx_relaxed,'ubx',ubx_relaxed,'lbg',lbg_relaxed,'ubg',ubg_relaxed);
                        % Cpu times
                        cpu_time_presolve_nlp_ii = toc(t_presolve_nlp_iter);
                        stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_presolve_nlp_ii ];
                        stats_nlp = solver_relaxed.stats();
                        nlp_iters_ii = stats_nlp.iter_count;
                        % Extract results and compute objective, infeasibility ect.
                        x_k = full(results_nlp.x);
                        lambda_x_k = full(results_nlp.lam_x);
                        h_comp_ii = full(mpec_casadi.h_comp_fun(x_k(1:dims.n_primal),p0));
                        h_std_ii= full(mpec_casadi.h_std_fun(x_k(1:dims.n_primal),p0));
                        f_opt_ii = full(mpec_casadi.f_fun(x_k(1:dims.n_primal),p0));
                        stats.iter.X_outer = [stats.iter.X_outer, x_k(1:dims.n_primal)];
                        if settings.verbose_solver
                            if strcmp(settings.relax_and_project_homotopy_parameter_steering,"Direct")
                                print_iter_stats(1,ii,f_opt_ii,h_std_ii,h_comp_ii,'NLP (Reg)',nlp_iters_ii,stats_nlp.return_status,p0_relaxed(end),norm(x_k_init-x_k),cpu_time_presolve_nlp_ii,1)
                            else
                                print_iter_stats(1,ii,f_opt_ii,h_std_ii,h_comp_ii,'NLP (Pen)',nlp_iters_ii,stats_nlp.return_status,p0_relaxed(end),norm(x_k_init-x_k),cpu_time_presolve_nlp_ii,1)
                            end
                        end
                        if isequal(stats_nlp.return_status,'Infeasible_Problem_Detected')
                            stats.phase_i_infeasibility_detected = true;
                            if settings.stop_if_nlp_infeasible
                                stats.problem_infeasible  = true;
                                stats.success = false;
                                stats.success_phase_i = false;
                                stats.stopping_criterion_fullfiled = true;
                                stats.solver_message = 'MPEC is locally infeasible.';
                            end
                        end
                        p0_relaxed(end) = settings.relax_and_project_kappa*p0_relaxed(end);
                        ii = ii+1;
                        x_k_init = x_k;
                        stats.n_nlp_total = stats.n_nlp_total+1;
                        x_k = x_k_init(1:dims.n_primal);
                    end
                end
                y_lpec_k_l = abs(x_k(dims.ind_x1))>=abs(x_k(dims.ind_x2)); % inital guess for active set/binaries
                tol_active_set = max(settings.tol_active,min(h_comp_ii,p0_relaxed(end)));

              case {'FeasibilityEllInfGeneral','FeasibilityEll1General'}
                n_g = length(mpec_casadi.g);
                % chekc if inital point already feasible
                h_total_k = full(mpec_casadi.h_total_fun(x_k,p0));
                if h_total_k <= settings.tol_feasibility
                    % TODO: add this maybe before any possible Phase I method;
                    stats.feasible_bnlp_found = true;
                    stats.success_phase_i = true;
                    x_k = project_to_bounds(x_k,lbx,ubx,dims);
                    fprintf('\n MPECopt: initial guess already feasible. \n')
                else
                    if n_g > 0
                        if strcmp(settings.initialization_strategy,"FeasibilityEll1General")
                            s = SX.sym('s',n_g); %define slack variables for ell_1 norm of generla constraints;
                            n_slacks = n_g;
                        else
                            s = SX.sym('s',1); %define slack variables for ell_inf norm of generla constraints;
                            n_slacks = 1;
                        end
                        % do a for loop for index updates to preserve the structure to some extent
                        g_s = [];
                        lbg_s = [];
                        ubg_s = [];
                        % (TODO: vectorize the for loop above for better perfomance)
                        for ii = 1:n_g
                            if strcmp(settings.initialization_strategy,"FeasibilityEll1General")
                                s_ii = s(ii);
                            else
                                s_ii = s;
                            end
                            % lower bound;
                            if solver_initialization.ubg(ii)~= -inf
                                g_s = [g_s; mpec_casadi.g(ii)+s_ii];
                                lbg_s = [lbg_s; solver_initialization.lbg(ii)];
                                ubg_s = [ubg_s; inf];
                            end
                            % upper bound;
                            if solver_initialization.ubg(ii)~= inf
                                g_s = [g_s; mpec_casadi.g(ii)-s_ii];
                                lbg_s = [lbg_s; -inf];
                                ubg_s = [ubg_s; solver_initialization.ubg(ii)];
                            end
                        end
                        if settings.feasibility_project_to_bounds
                            x_projected = project_to_bounds(solver_initialization.x0,solver_initialization.lbx,solver_initialization.ubx,dims);
                            solver_initialization.x0 = x_projected;
                        end
                        s_max = full(mpec_casadi.h_std_fun(solver_initialization.x0,solver_initialization.p0))+1e-1;
                        s_lb = settings.feasibility_s_lower_bound*ones(n_slacks,1);
                        s_ub = 1.1*s_max*ones(n_slacks,1);
                        s0 = 0.9*s_max*ones(n_slacks,1);

                        % new mpec for feasibility
                        mpec_feas = struct;
                        mpec_feas.f = sum(s);
                        mpec_feas.g = g_s;
                        mpec_feas.x = [s;mpec_casadi.x];
                        mpec_feas.G = mpec_casadi.x1;
                        mpec_feas.H = mpec_casadi.x2;
                        solver_initialization_feas = solver_initialization;
                        solver_initialization_feas.lbg  = lbg_s;
                        solver_initialization_feas.ubg = ubg_s;
                        solver_initialization_feas.lbx = [s_lb;solver_initialization.lbx];
                        solver_initialization_feas.ubx  = [s_ub; solver_initialization.ubx];
                        solver_initialization_feas.x0 = [s0;solver_initialization.x0];
                        [mpec_feas_casadi, dims_feas, solver_initialization_feas] =  create_mpec_functions(mpec_feas,solver_initialization_feas,settings);
                        dims_feas.n_slacks = n_slacks;
                        lpec_feas_casadi = create_lpec_functions(mpec_feas_casadi,dims_feas,settings,solver_initialization_feas);
                        solver_initialization_feas.y_lpec_k_l = y_lpec_k_l;
                        solver_initialization_feas.d_lpec_k_l = nan;
                        phase_ii = false;
                        [solution_feas,stats_feas] = mpecopt_phase_ii(mpec_feas_casadi,lpec_feas_casadi,dims_feas,settings,solver_initialization_feas,stats,phase_ii);
                        s_k = solution_feas.x(1:n_slacks);
                        % get cpu times from phase i
                        stats.iter.cpu_time_lpec_phase_i_iter = stats_feas.iter.cpu_time_lpec_phase_i_iter;
                        stats.iter.cpu_time_nlp_phase_i_iter = stats_feas.iter.cpu_time_nlp_phase_i_iter;
                        s_infeasiblity = max(s_k);
                        if s_infeasiblity >= settings.tol_feasibility
                            stats.stopping_criterion_fullfiled = true;
                            stats.problem_infeasible  = true;
                            stats.success = false;
                            stats.success_phase_i = false;
                            stats.solver_message = 'MPEC is locally infeasible.';
                            x_k = solution_feas.x(1+n_slacks:end);
                        else
                            stats.feasible_bnlp_found = true;
                            stats.success_phase_i = true;
                            x_k = solution_feas.x(n_slacks+1:end);
                            if full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)) > 1e-3 && settings.debug_mode_on
                                keyboard;
                            end
                        end
                    else
                        stats.feasible_bnlp_found = true;
                        stats.success_phase_i = true;
                        x_k = project_to_bounds(x_k,lbx,ubx,dims);
                        fprintf('\n MPECopt: MPEC has only complementarity and bound constraints, feasible point found by projection to bounds.\n')
                    end
                end

              case 'RandomActiveSet' % just take the user provided x0;
                                     % rng(1,"twister"); % to have reproducable problem sets
                ii = 1;
                while (ii <= settings.max_recovery_iters) && ~stats.problem_infeasible
                    I_0_plus = boolean(round(rand(dims.n_comp,1)));
                    I_plus_0 = ~I_0_plus;
                    lbx_bnlp_k = lbx;
                    ubx_bnlp_k = ubx;
                    ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                    ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                    % solve BNLP/NLP
                    x_k_init = x_k;
                    t_presolve_nlp_iter = tic;
                    results_nlp = mpec_casadi.solver('x0',x_k,'p',p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',lbg,'ubg',ubg);
                    cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                    stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k ];
                    x_k = full(results_nlp.x); lambda_x_k  = full(results_nlp.lam_x);
                    stats_nlp = mpec_casadi.solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                    h_comp_k = full(mpec_casadi.h_comp_fun(x_k,p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,p0)); f_opt_k = full(results_nlp.f);
                    if settings.verbose_solver
                        print_iter_stats(1,ii,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                    end
                    ii = ii+1;
                    stats.iter.X_outer = [stats.iter.X_outer, x_k];
                    if ismember(stats_nlp.return_status, {'Solve_Succeeded', 'Search_Direction_Becomes_Too_Small', 'Solved_To_Acceptable_Level'})
                        stats.feasible_bnlp_found = true;
                        stats.success_phase_i = true;
                        break;
                    end
                end
                % TODO: IF INFEASIBLE QUIT
              case 'AllBiactive'
                lbx_bnlp_k = lbx; ubx_bnlp_k = ubx; x_k_init = x_k;
                lbx_bnlp_k(dims.ind_x1) = 0; ubx_bnlp_k(dims.ind_x1) = 0; lbx_bnlp_k(dims.ind_x2) = 0; ubx_bnlp_k(dims.ind_x2) = 0;
                t_presolve_nlp_iter = tic;
                results_nlp = mpec_casadi.solver('x0', x_k, 'p',p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',lbg,'ubg',ubg); % solve BNLP/NLP
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k];
                x_k = full(results_nlp.x);
                lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = mpec_casadi.solver.stats();
                nlp_iters_k = stats_nlp.iter_count;
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k,p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,p0)); f_opt_k = full(results_nlp.f);
                if settings.verbose_solver
                    print_iter_stats('BNLP',1,f_opt_k,h_std_k,h_comp_k,'NLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.n_nlp_total = stats.n_nlp_total + 1;

                % TODO: IF INFEASIBLE QUIT
              case 'TakeInitialGuessDirectly'
                if settings.project_guess_to_bounds
                    x_k = project_to_bounds(x_k,lbx,ubx,dims);
                end
              case 'TakeInitialGuessActiveSet'
                % Make active set guess from x_k and solve corresponding BNLP;
                I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                lbx_bnlp_k = lbx;
                ubx_bnlp_k = ubx;
                ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                % solve BNLP/NLP
                x_k_init = x_k;
                t_presolve_nlp_iter = tic;
                results_nlp = mpec_casadi.solver('x0',x_k,'p',p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',lbg,'ubg',ubg);
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k ];
                x_k = full(results_nlp.x);
                lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = mpec_casadi.solver.stats();
                nlp_iters_k = stats_nlp.iter_count;
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k,p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,p0)); f_opt_k = full(results_nlp.f);
                if settings.verbose_solver
                    print_iter_stats(1,1,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.n_nlp_total = stats.n_nlp_total + 1;
                % stats.succes+p

              case 'TakeProvidedActiveSet'
                if isfield('y0',solver_initialization)
                    if isequal(length(y0),dims.n_comp)
                        y_lpec_k_l = solver_initialization.y0;
                        I_plus_0 = y_lpec_k_l ==1;
                        I_0_plus = y_lpec_k_l == 0;
                    else
                        warning('\n TakeProvidedActiveSet: lenth of y0 is %d, but required is %d',length(solver_initialization.y0))
                        I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                        I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                    end
                else
                    warning('\n TakeProvidedActiveSet as initialization strategy chosen, but no y0 provided. Resporting to TakeInitialGuessActiveSet')
                    I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                    I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                end
                lbx_bnlp_k = lbx;
                ubx_bnlp_k = ubx;
                ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                % solve BNLP/NLP
                x_k_init = x_k;
                t_presolve_nlp_iter = tic;
                results_nlp = solver('x0',x_k,'p',p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',lbg,'ubg',ubg);
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter; cpu_time_bnlp_k];
                x_k = full(results_nlp.x); lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                h_comp_k = full(h_comp_fun(x_k,p0)); h_std_k = full(h_std_fun(x_k,p0)); f_opt_k = full(results_nlp.f);
                if settings.verbose_solver
                    print_iter_stats(1,1,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.success_phase_i = true;
            end
            stats.rho_TR_final = rho_TR_k_l;
            %  ---- make sure feasible point is declared sucessful ---
            if full(mpec_casadi.h_total_fun(x_k,p0)) <= settings.tol_feasibility
                stats.feasible_bnlp_found = true; stats.success_phase_i = true;
            else
                stats.feasible_bnlp_found = false; stats.success_phase_i = false;
            end
            stats.phase_i_infeasibility_detected = ~stats.success_phase_i;

            % -------------------- Check is the BNLP solution S-stationary ------------------------------
            h_total_phase_i  = full(mpec_casadi.h_total_fun(x_k,p0));   % infeasibility

            if stats.success_phase_i && settings.stop_if_S_stationary && ~stats.stopping_criterion_fullfiled && h_total_phase_i <= settings.tol_feasibility && ~(strcmp(settings.initialization_strategy,'FeasibilityEll1General') || strcmp(settings.initialization_strategy,'FeasibilityEllInfGeneral'))
                active_set_estimate_k = find_active_sets(x_k, dims, settings.tol_active); % check is it S-stationary;
                if sum(active_set_estimate_k.I_00) == 0
                    stats.stopping_criterion_fullfiled = true;
                    stats.multiplier_based_stationarity = 'S';
                    stats.solver_message = 'S-stationary point found sucessfully in Phase I.';
                    stats.solved_in_phase_i = true;
                    stats.b_stationarity = true;
                    stats.success = true;
                    stats.f_lpec = 0;
                    % todo: resolve a tnlp to check?
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
                        [multiplier_based_stationarity, solver_message] = determine_multipliers_based_stationary_point(x_k,lambda_x_k,dims,settings);
                        if strcmp(stats.multiplier_based_stationarity,'S')
                            stats.stopping_criterion_fullfiled = true;
                            stats.solver_message = 'S-stationary point found sucessfully in presolve.';
                            stats.success = true;
                            stats.b_stationarity = true;
                            stats.solved_in_phase_i = true;
                            stats.f_lpec = 0;
                        end
                    end
                end
            end
            stats.cpu_time_phase_i = toc(t_phase_i_start);
            % ------------------------------- end of Phase I -----------------------------------------------

            %% ---------------------- Main optimization loop of Phase II -----------------------------------
            solver_initialization.x0 = x_k;
            solver_initialization.y_lpec_k_l = y_lpec_k_l;
            solver_initialization.d_lpec_k_l = d_lpec_k_l;
            phase_ii = true;
            if settings.verbose_solver
                print_phase_ii();
                print_iter_header();
            end
            [solution,stats] = mpecopt_phase_ii(mpec_casadi,lpec_casadi,dims,settings,solver_initalization,stats,phase_ii);

            if  stats.solved_in_phase_i
                stats.cpu_time_phase_ii = 1e-10;
            end

            %% ------------------------ Final verbose ----------------------------------
            fprintf('\n');
            if settings.verbose_summary
                print_iter_summary(solution.f,stats.h_std,stats.comp_res,stats.solver_message,stats.multiplier_based_stationarity,stats.b_stationarity,stats.n_biactive,stats.f_lpec,stats.rho_TR_final);
                fprintf('\n');
                if settings.verbose_extended_summary
                    print_iter_details(stats.n_nlp_total,stats.n_lpec_total,stats.iter.rho_TR_iter(end),cpu_time_prepare_mpec+cpu_time_generate_nlp_solvers,stats.cpu_time_total,stats.cpu_time_phase_i,stats.cpu_time_main_loop)
                    try
                        fprintf('LPEC objective...............:\t %2.2e\n',results_lpec.f_opt);
                    catch
                    end
                    fprintf('Final TR radius..............:\t %2.2e\n',rho_TR_k_l);
                    fprintf('Equality constraints.........:\t %d\n',stats.n_eq);
                    fprintf('Simple eq. constraints.......:\t %d\n',stats.n_box_simple);
                    fprintf('Active inequalites...........:\t %d/%d\n',stats.n_active_ineq,stats.n_ineq);
                    fprintf('Active box constraints.......:\t %d/%d\n',stats.n_active_box,stats.n_box);
                end
            end


            %% -----------------------------Detalied solver stats---------------------------
            % summarized more specific
            try
                stats.cpu_time_total = stats.cpu_time_phase_i+stats.cpu_time_phase_ii;
            catch
                stats.cpu_time_total = stats.cpu_time_phase_ii;
            end


            stats.cpu_time_prepare_mpec = cpu_time_prepare_mpec;
            stats.cpu_time_lpec_preparation = sum(stats.iter.cpu_time_lpec_preparation_iter);
            stats.cpu_time_lpec_phase_i = sum(stats.iter.cpu_time_lpec_phase_i_iter);
            stats.cpu_time_lpec_phase_ii = sum(stats.iter.cpu_time_lpec_phase_ii_iter);
            stats.cpu_time_lpec = stats.cpu_time_lpec_phase_i+stats.cpu_time_lpec_phase_ii;
            stats.cpu_time_nlp_phase_i = sum(stats.iter.cpu_time_nlp_phase_i_iter);
            stats.cpu_time_nlp_phase_ii = sum(stats.iter.cpu_time_nlp_phase_ii_iter);
            stats.cpu_time_nlp = stats.cpu_time_nlp_phase_i+stats.cpu_time_nlp_phase_ii;
            stats.cpu_time_globalization = sum(stats.iter.cpu_time_globalization_iter);
            if isempty(stats.iter.cpu_time_nlp_iter)
                stats.iter.cpu_time_nlp_iter = nan;
            end
            stats.dims = dims;

            %% update info
            obj.stats = stats;
        end
    end
end
