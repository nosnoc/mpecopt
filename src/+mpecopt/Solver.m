classdef Solver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        opts % all options of mpecopt
        mpec % casadi expressions defining the mpec/mpcc
        mpec_casadi % all casadi functions and solvers needed for the mpec/mpec
        lpec_casadi % all casadi functions and solvers needed for the lpec/lpec
        dims % all problem dimensions and relevant index sets
        solver_initialization % solver initialization data, x0, lower and upper bounds of variables and constraints
        stats % solution statistics, mostly cpu times, and some qualitative solution information
    end

    methods
        function obj = Solver(mpec, opts)
            t_prepare_mpec = tic;
            obj.opts = opts;
            obj.mpec = mpec;
            obj.create_mpec_functions();
            obj.solver_initialization = struct();
            obj.create_lpec_functions();
            
            cpu_time_prepare_mpec = toc(t_prepare_mpec);
            obj.stats.cpu_time_prepare_mpec = cpu_time_prepare_mpec;
        end

        function [solution,stats] = solve(obj, solver_initialization)
            arguments
                obj
                solver_initialization = []
            end
            %% Get data
            import casadi.*
            obj.process_solver_initialization(solver_initialization);
            solver_initialization = obj.solver_initialization;
            mpec = obj.mpec;
            opts = obj.opts;
            mpec = obj.mpec;
            mpec_casadi = obj.mpec_casadi;
            lpec_casadi = obj.lpec_casadi;
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


            rho_TR_k_l = opts.rho_TR_phase_i_init;
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

            if opts.plot_mpec_multipliers
                figure;     % TODO: add figure handle that is passed to the multipl function.
            end

            % ------------------ Prepare homotopy solver for Phase I -----------------------------
            if strcmp(opts.initialization_strategy,"RelaxAndProject")
                t_generate_nlp_solvers = tic;
                % TODO : input just mpec_casadi and solver_initialization
                [solver_relaxed,x_k_init,p0_relaxed,lbx_relaxed,ubx_relaxed,lbg_relaxed,ubg_relaxed] = create_phase_i_nlp_solver(mpec_casadi.f,mpec_casadi.g,mpec_casadi.x,mpec_casadi.x1,mpec_casadi.x2,mpec_casadi.p,solver_initialization.lbx,solver_initialization.ubx,solver_initialization.lbg,solver_initialization.ubg,solver_initialization.p0,x_k,opts,dims);
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
            % solver = nlpsol('solver', 'ipopt', nlp, opts.settings_casadi_nlp);
            % stats.cpu_time_generate_nlp_solvers = toc(t_generate_nlp_solvers);
            % mpec_casadi.solver = solver;

            %%  ---------------------- Phase I - compute feasible point ------------------------
            if opts.verbose_solver
                print_phase_i();
                print_iter_header();
            end
            h_comp_ii = full(mpec_casadi.h_comp_fun(x_k(1:dims.n_primal),solver_initialization.p0));
            h_std_ii = full(mpec_casadi.h_std_fun(x_k(1:dims.n_primal),solver_initialization.p0));
            f_opt_ii = full(mpec_casadi.f_fun(x_k(1:dims.n_primal),solver_initialization.p0));
            if opts.verbose_solver
                print_iter_stats('I',0,f_opt_ii,h_std_ii,h_comp_ii,'/',0,'Initial guess',0,0,0,1)
            end

            t_phase_i_start = tic;
            switch opts.initialization_strategy
              case 'RelaxAndProject'
                h_comp_ii = 1;
                ii = 1;
                y_lpec_k_l = double(x_k(dims.ind_x1)> x_k(dims.ind_x2));
                d_lpec_k_l = zeros(dims.n_primal,1);
                while (ii <= opts.max_recovery_iters) && ~stats.problem_infeasible && ~feasible_bnlp_found
                    switch opts.bnlp_projection_strategy
                      case {'LPEC'}
                        if opts.lpec_solve_if_comp_feasible
                            rho_min = 1.01*max(abs(min(x_k(dims.ind_x1),x_k(dims.ind_x2))));
                            if rho_min > opts.rho_TR_phase_i_init || ii == 1
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
                            lpec = create_lpec_subproblem(x_k,solver_initialization.p0,rho_TR_k_l,lpec_casadi,dims,opts,opts.tol_active);
                            stats.iter.cpu_time_lpec_preparation_iter = [stats.iter.cpu_time_lpec_preparation_iter;toc(t_lpec_preparation_iter)];
                            %  Initial guess and TR for the LPEC
                            y_lpec_k_previous = y_lpec_k_l; % to keep track of active set chnages
                            if opts.warm_start_lpec_phase_i
                                lpec.d_lpec = d_lpec_k_l;
                                lpec.y_lpec = y_lpec_k_l; % inital guess for bin. variablels.
                            end

                            rho_TR_k_l_lb = opts.realx_and_project_scale_factor_rho_tr*max(abs(min(x_k(dims.ind_x1),x_k(dims.ind_x2))));
                            rho_TR_k_l_ub = opts.realx_and_project_scale_factor_rho_tr*max(max(abs(x_k(dims.ind_x1)),abs(x_k(dims.ind_x2))));

                            if opts.relax_and_project_tighter_TR
                                rho_TR_k_l = min(rho_TR_k_l_lb, opts.rho_TR_phase_i_init);
                                rho_TR_k_l = max(1e-4,rho_TR_k_l); % to avoid unerasonably small TR at almost fesaible points;
                                                                   % [rho_TR_k_l_lb,  p0_relaxed(end) ]
                            else
                                rho_TR_k_l = max(rho_TR_k_l_lb, opts.rho_TR_phase_i_init);
                            end
                            lpec.rho_TR  =  rho_TR_k_l;
                            stats.iter.rho_TR_iter = [stats.iter.rho_TR_iter, rho_TR_k_l]; % store TR radius
                                                                                           % Solve LPEC
                            [results_lpec,stats_lpec] = lpec_solver(lpec,opts.settings_lpec);
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
                            h_comp_lpec_k_l = full(mpec_casadi.h_comp_fun(x_trail_lpec,solver_initialization.p0));
                            h_std_lpec_k_l = full(mpec_casadi.h_std_fun(x_trail_lpec,solver_initialization.p0));
                            if opts.verbose_solver
                                print_iter_stats(1,ii,f_lin_opt_k_l,h_std_lpec_k_l,h_comp_lpec_k_l,'LPEC',stats_lpec.nodecount,stats_lpec.solver_message,lpec.rho_TR,norm(d_lpec_k_l),stats_lpec.cpu_time,' ')
                            end
                            if opts.use_lpec_fallback_strategy_phase_i  && ii > ceil(opts.lpec_recovery_start*opts.max_recovery_iters)
                                if strcmp(stats_lpec.solver_message,'INFEASIBLE') || (strcmp(stats_lpec.solver_message,'NODE_LIMIT') && isnan(results_lpec.f_opt))
                                    if opts.debug_mode_on
                                        keyboard;
                                    end
                                    [results_lpec,stats_lpec] = lpec_fallback_strategy(lpec,opts,results_lpec,stats_lpec);
                                    stats.n_lpec_total = stats.n_lpec_total + results_lpec.n_lpecs_solved;
                                    lpec_solution_exists = stats_lpec.lpec_solution_exists;
                                    d_lpec_k_l = results_lpec.d_lpec;  y_lpec_k_l = results_lpec.y_lpec; f_lin_opt_k_l = results_lpec.f_opt;
                                    x_trail_lpec = x_k + d_lpec_k_l;
                                end
                            end
                            %
                            if lpec_solution_exists
                                % --------------------------- Check if B-stationary point found --------------------------
                                h_total_k = full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0));
                                nabla_f_k = full(mpec_casadi.nabla_f_fun(x_k,solver_initialization.p0));
                                if (h_total_k <= opts.tol) && ((abs(f_lin_opt_k_l) <= opts.tol_B_stationarity || norm(nabla_f_k) <= opts.tol_B_stationarity))  % if objective zero (either if cost gradient zero, or solution leads to it) = then set step to zero => B stationarity
                                    if opts.reset_lpec_objective
                                        d_lpec_k_l = d_lpec_k_l*0; % if the current point is feasible, and the objective is zero, then d = 0 is also a solution of the lpec (occurs if a solution is not on the verties of the lp)
                                        f_lin_opt_k_l = 0;
                                    end
                                end
                                if norm(d_lpec_k_l) <= opts.tol_B_stationarity
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
                    h_total_k = full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0));
                    if h_total_k <= opts.tol_feasibility && ii>1
                        feasible_bnlp_found = true;
                    end

                    if active_set_guess_exists && ~feasible_bnlp_found
                        lbx_bnlp_k = solver_initialization.lbx;
                        ubx_bnlp_k = solver_initialization.ubx;
                        ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                        ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                        t_presolve_nlp_iter = tic;
                        results_nlp = mpec_casadi.solver('x0',x_k,'p',solver_initialization.p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',solver_initialization.lbg,'ubg',solver_initialization.ubg);
                        cpu_time_bnlp_ii = toc(t_presolve_nlp_iter);
                        stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter; cpu_time_bnlp_ii];
                        x_k = full(results_nlp.x);
                        lambda_x_k  = full(results_nlp.lam_x);
                        stats_nlp = mpec_casadi.solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                        h_comp_k = full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0));
                        h_std_k = full(mpec_casadi.h_std_fun(x_k,solver_initialization.p0));
                        f_opt_k = full(results_nlp.f);
                        if opts.verbose_solver
                            print_iter_stats(1,ii,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,rho_TR_k_l,norm(x_k_init(1:dims.n_primal)-x_k),cpu_time_bnlp_ii ,1)
                        end
                        if full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0)) <= opts.tol_feasibility
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
                        h_comp_ii = full(mpec_casadi.h_comp_fun(x_k(1:dims.n_primal),solver_initialization.p0));
                        h_std_ii= full(mpec_casadi.h_std_fun(x_k(1:dims.n_primal),solver_initialization.p0));
                        f_opt_ii = full(mpec_casadi.f_fun(x_k(1:dims.n_primal),solver_initialization.p0));
                        stats.iter.X_outer = [stats.iter.X_outer, x_k(1:dims.n_primal)];
                        if opts.verbose_solver
                            if strcmp(opts.relax_and_project_homotopy_parameter_steering,"Direct")
                                print_iter_stats(1,ii,f_opt_ii,h_std_ii,h_comp_ii,'NLP (Reg)',nlp_iters_ii,stats_nlp.return_status,p0_relaxed(end),norm(x_k_init-x_k),cpu_time_presolve_nlp_ii,1)
                            else
                                print_iter_stats(1,ii,f_opt_ii,h_std_ii,h_comp_ii,'NLP (Pen)',nlp_iters_ii,stats_nlp.return_status,p0_relaxed(end),norm(x_k_init-x_k),cpu_time_presolve_nlp_ii,1)
                            end
                        end
                        if isequal(stats_nlp.return_status,'Infeasible_Problem_Detected')
                            stats.phase_i_infeasibility_detected = true;
                            if opts.stop_if_nlp_infeasible
                                stats.problem_infeasible  = true;
                                stats.success = false;
                                stats.success_phase_i = false;
                                stats.stopping_criterion_fullfiled = true;
                                stats.solver_message = 'MPEC is locally infeasible.';
                            end
                        end
                        p0_relaxed(end) = opts.relax_and_project_kappa*p0_relaxed(end);
                        ii = ii+1;
                        x_k_init = x_k;
                        stats.n_nlp_total = stats.n_nlp_total+1;
                        x_k = x_k_init(1:dims.n_primal);
                    end
                end
                y_lpec_k_l = abs(x_k(dims.ind_x1))>=abs(x_k(dims.ind_x2)); % inital guess for active set/binaries
                tol_active_set = max(opts.tol_active,min(h_comp_ii,p0_relaxed(end)));

              case {'FeasibilityEllInfGeneral','FeasibilityEll1General'}
                n_g = length(mpec_casadi.g);
                % chekc if inital point already feasible
                h_total_k = full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0));
                if h_total_k <= opts.tol_feasibility
                    % TODO: add this maybe before any possible Phase I method;
                    stats.feasible_bnlp_found = true;
                    stats.success_phase_i = true;
                    x_k = project_to_bounds(x_k,solver_initialization.lbx,solver_initialization.ubx,dims);
                    fprintf('\n MPECopt: initial guess already feasible. \n')
                else
                    if n_g > 0
                        if strcmp(opts.initialization_strategy,"FeasibilityEll1General")
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
                            if strcmp(opts.initialization_strategy,"FeasibilityEll1General")
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
                        if opts.feasibility_project_to_bounds
                            x_projected = project_to_bounds(solver_initialization.x0,solver_initialization.lbx,solver_initialization.ubx,dims);
                            solver_initialization.x0 = x_projected;
                        end
                        s_max = full(mpec_casadi.h_std_fun(solver_initialization.x0,solver_initialization.p0))+1e-1;
                        s_lb = opts.feasibility_s_lower_bound*ones(n_slacks,1);
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
                        [mpec_feas_casadi, dims_feas, solver_initialization_feas] =  create_mpec_functions(mpec_feas,solver_initialization_feas,opts);
                        dims_feas.n_slacks = n_slacks;
                        lpec_feas_casadi = create_lpec_functions(mpec_feas_casadi,dims_feas,opts,solver_initialization_feas);
                        solver_initialization_feas.y_lpec_k_l = y_lpec_k_l;
                        solver_initialization_feas.d_lpec_k_l = nan;
                        phase_ii = false;
                        [solution_feas,stats_feas] = mpecopt_phase_ii(mpec_feas_casadi,lpec_feas_casadi,dims_feas,opts,solver_initialization_feas,stats,phase_ii);
                        s_k = solution_feas.x(1:n_slacks);
                        % get cpu times from phase i
                        stats.iter.cpu_time_lpec_phase_i_iter = stats_feas.iter.cpu_time_lpec_phase_i_iter;
                        stats.iter.cpu_time_nlp_phase_i_iter = stats_feas.iter.cpu_time_nlp_phase_i_iter;
                        s_infeasiblity = max(s_k);
                        if s_infeasiblity >= opts.tol_feasibility
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
                            if full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)) > 1e-3 && opts.debug_mode_on
                                keyboard;
                            end
                        end
                    else
                        stats.feasible_bnlp_found = true;
                        stats.success_phase_i = true;
                        x_k = project_to_bounds(x_k,solver_initialization.lbx,solver_initialization.ubx,dims);
                        fprintf('\n MPECopt: MPEC has only complementarity and bound constraints, feasible point found by projection to bounds.\n')
                    end
                end

              case 'RandomActiveSet' % just take the user provided x0;
                                     % rng(1,"twister"); % to have reproducable problem sets
                ii = 1;
                while (ii <= opts.max_recovery_iters) && ~stats.problem_infeasible
                    I_0_plus = boolean(round(rand(dims.n_comp,1)));
                    I_plus_0 = ~I_0_plus;
                    lbx_bnlp_k = solver_initialization.lbx;
                    ubx_bnlp_k = solver_initialization.ubx;
                    ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                    ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                    % solve BNLP/NLP
                    x_k_init = x_k;
                    t_presolve_nlp_iter = tic;
                    results_nlp = mpec_casadi.solver('x0',x_k,'p',solver_initialization.p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',solver_initialization.lbg,'ubg',solver_initialization.ubg);
                    cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                    stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k ];
                    x_k = full(results_nlp.x); lambda_x_k  = full(results_nlp.lam_x);
                    stats_nlp = mpec_casadi.solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                    h_comp_k = full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,solver_initialization.p0)); f_opt_k = full(results_nlp.f);
                    if opts.verbose_solver
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
                lbx_bnlp_k = solver_initialization.lbx; ubx_bnlp_k = solver_initialization.ubx; x_k_init = x_k;
                lbx_bnlp_k(dims.ind_x1) = 0; ubx_bnlp_k(dims.ind_x1) = 0; lbx_bnlp_k(dims.ind_x2) = 0; ubx_bnlp_k(dims.ind_x2) = 0;
                t_presolve_nlp_iter = tic;
                results_nlp = mpec_casadi.solver('x0', x_k, 'p',solver_initialization.p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',solver_initialization.lbg,'ubg',solver_initialization.ubg); % solve BNLP/NLP
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k];
                x_k = full(results_nlp.x);
                lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = mpec_casadi.solver.stats();
                nlp_iters_k = stats_nlp.iter_count;
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,solver_initialization.p0)); f_opt_k = full(results_nlp.f);
                if opts.verbose_solver
                    print_iter_stats('BNLP',1,f_opt_k,h_std_k,h_comp_k,'NLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.n_nlp_total = stats.n_nlp_total + 1;

                % TODO: IF INFEASIBLE QUIT
              case 'TakeInitialGuessDirectly'
                if opts.project_guess_to_bounds
                    x_k = project_to_bounds(x_k,solver_initialization.lbx,solver_initialization.ubx,dims);
                end
              case 'TakeInitialGuessActiveSet'
                % Make active set guess from x_k and solve corresponding BNLP;
                I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                lbx_bnlp_k = solver_initialization.lbx;
                ubx_bnlp_k = solver_initialization.ubx;
                ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                % solve BNLP/NLP
                x_k_init = x_k;
                t_presolve_nlp_iter = tic;
                results_nlp = mpec_casadi.solver('x0',x_k,'p',solver_initialization.p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',solver_initialization.lbg,'ubg',solver_initialization.ubg);
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter;cpu_time_bnlp_k ];
                x_k = full(results_nlp.x);
                lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = mpec_casadi.solver.stats();
                nlp_iters_k = stats_nlp.iter_count;
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,solver_initialization.p0)); f_opt_k = full(results_nlp.f);
                if opts.verbose_solver
                    print_iter_stats(1,1,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.n_nlp_total = stats.n_nlp_total + 1;
                % stats.succes+p

              case 'TakeProvidedActiveSet'
                if isfield(solver_initialization,'y0')
                    if isequal(length(solver_initialization.y0),dims.n_comp)
                        y_lpec_k_l = solver_initialization.y0;
                        I_plus_0 = y_lpec_k_l ==1;
                        I_0_plus = y_lpec_k_l == 0;
                    else
                        warning('TakeProvidedActiveSet: lenth of y0 is %d, but required is %d',length(solver_initialization.y0))
                        I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                        I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                    end
                else
                    warning('\n TakeProvidedActiveSet as initialization strategy chosen, but no y0 provided. Resporting to TakeInitialGuessActiveSet')
                    I_0_plus = x_k(dims.ind_x1)<=x_k(dims.ind_x2);
                    I_plus_0 = x_k(dims.ind_x1)>x_k(dims.ind_x2);
                end
                lbx_bnlp_k = solver_initialization.lbx;
                ubx_bnlp_k = solver_initialization.ubx;
                ubx_bnlp_k(dims.ind_x1(I_0_plus)) = 0;
                ubx_bnlp_k(dims.ind_x2(I_plus_0)) = 0;
                % solve BNLP/NLP
                x_k_init = x_k;
                t_presolve_nlp_iter = tic;
                results_nlp = mpec_casadi.solver('x0',x_k,'p',solver_initialization.p0,'lbx',lbx_bnlp_k,'ubx',ubx_bnlp_k,'lbg',solver_initialization.lbg,'ubg',solver_initialization.ubg);
                cpu_time_bnlp_k = toc(t_presolve_nlp_iter);
                stats.iter.cpu_time_nlp_phase_i_iter = [stats.iter.cpu_time_nlp_phase_i_iter; cpu_time_bnlp_k];
                x_k = full(results_nlp.x); lambda_x_k  = full(results_nlp.lam_x);
                stats_nlp = mpec_casadi.solver.stats(); nlp_iters_k = stats_nlp.iter_count; stats.n_nlp_total = stats.n_nlp_total + 1;
                h_comp_k = full(mpec_casadi.h_comp_fun(x_k,solver_initialization.p0)); h_std_k = full(mpec_casadi.h_std_fun(x_k,solver_initialization.p0)); f_opt_k = full(results_nlp.f);
                if opts.verbose_solver
                    print_iter_stats(1,1,f_opt_k,h_std_k,h_comp_k,'BNLP',nlp_iters_k,stats_nlp.return_status,nan,norm(x_k_init-x_k),cpu_time_bnlp_k,1)
                end
                stats.iter.X_outer = [stats.iter.X_outer, x_k];
                stats.success_phase_i = true;
            end
            stats.rho_TR_final = rho_TR_k_l;
            %  ---- make sure feasible point is declared sucessful ---
            if full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0)) <= opts.tol_feasibility
                stats.feasible_bnlp_found = true; stats.success_phase_i = true;
            else
                stats.feasible_bnlp_found = false; stats.success_phase_i = false;
            end
            stats.phase_i_infeasibility_detected = ~stats.success_phase_i;

            % -------------------- Check is the BNLP solution S-stationary ------------------------------
            h_total_phase_i  = full(mpec_casadi.h_total_fun(x_k,solver_initialization.p0));   % infeasibility

            if stats.success_phase_i && opts.stop_if_S_stationary && ~stats.stopping_criterion_fullfiled && h_total_phase_i <= opts.tol_feasibility && ~(strcmp(opts.initialization_strategy,'FeasibilityEll1General') || strcmp(opts.initialization_strategy,'FeasibilityEllInfGeneral'))
                active_set_estimate_k = find_active_sets(x_k, dims, opts.tol_active); % check is it S-stationary;
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
                    if opts.resolve_tnlp_for_stationarity && ~opts.PieceNLPStartegy == "TNLP"
                        % get rid of bilinear constraints
                        ubg_relaxed(ind_comp) = +inf;
                        if strcmp(opts.relax_and_project_homotopy_parameter_steering,'Ell_1') || strcmp(opts.relax_and_project_homotopy_parameter_steering,'Ell_inf')
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
                        [multiplier_based_stationarity, solver_message] = determine_multipliers_based_stationary_point(x_k,lambda_x_k,dims,opts);
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
            if opts.verbose_solver
                print_phase_ii();
                print_iter_header();
            end
            [solution,stats] = mpecopt_phase_ii(mpec_casadi,lpec_casadi,dims,opts,solver_initialization,stats,phase_ii);

            if  stats.solved_in_phase_i
                stats.cpu_time_phase_ii = 1e-10;
            end

            %% ------------------------ Final verbose ----------------------------------
            fprintf('\n');
            if opts.verbose_summary
                print_iter_summary(solution.f,stats.h_std,stats.comp_res,stats.solver_message,stats.multiplier_based_stationarity,stats.b_stationarity,stats.n_biactive,stats.f_lpec,stats.rho_TR_final);
                fprintf('\n');
                if opts.verbose_extended_summary
                    print_iter_details(stats.n_nlp_total,stats.n_lpec_total,stats.iter.rho_TR_iter(end),stats.cpu_time_prepare_mpec+cpu_time_generate_nlp_solvers,stats.cpu_time_total,stats.cpu_time_phase_i,stats.cpu_time_main_loop)
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

        function out = cat(dim,varargin)
            nosnoc.error('invalid', 'Invalid Operation')
        end

        function varargout = size(obj,varargin)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            varargout = {1};
        end
        
        function ind = end(obj,k,n)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            ind = 1;
        end
    end

    methods(Access=protected)
        function varargout = parenReference(obj, index_op)
            % TODO(@anton) this returns mpccresults struct from homotopy
            import casadi.*;
            p = inputParser;
            addParameter(p, 'x0', []);
            addParameter(p, 'y0', []);
            addParameter(p, 'lbx', []);
            addParameter(p, 'ubx', []);
            addParameter(p, 'lbg', []);
            addParameter(p, 'ubg', []);
            addParameter(p, 'p', []);
            addParameter(p, 'lam_g0', []);
            addParameter(p, 'lam_x0', []);
            parse(p, index_op(1).Indices{:});

            dims = obj.dims;

            if ~isempty(p.Results.x0)
                obj.solver_initialization.x0 = p.Results.x0;
            else
                obj.solver_initialization.x0 = zeros(dims.n_x,1);
            end
            if ~isempty(p.Results.lbx)
                obj.solver_initialization.lbx = p.Results.lbx;
            else
                obj.solver_initialization.lbx = -inf(dims.n_x,1);
            end
            if ~isempty(p.Results.ubx)
                obj.solver_initialization.ubx = p.Results.ubx;
            else
                obj.solver_initialization.ubx = inf(dims.n_x,1);
            end
            if ~isempty(p.Results.lbg)
                obj.solver_initialization.lbg = p.Results.lbg;
            else
                obj.solver_initialization.lbg = zeros(dims.n_g,1);
            end
            if ~isempty(p.Results.ubg)
                obj.solver_initialization.ubg = p.Results.ubg;
            else
                obj.solver_initialization.ubg = zeros(dims.n_g,1);
            end
            if ~isempty(p.Results.p)
                obj.solver_initialization.p0 = p.Results.p;
            else
                obj.solver_initialization.p0 = zeros(dims.n_p,1);
            end
            if ~isempty(p.Results.lam_g0)
                obj.solver_initialization.lam_g0 = p.Results.lam_g0;
            else
                obj.solver_initialization.lam_g0 = zeros(dims.n_g, 1);
            end
            if ~isempty(p.Results.lam_x0)
                obj.solver_initialization.lam_x0 = p.Results.lam_x0;
            else
                obj.solver_initialization.lam_x0 = zeros(dims.n_x, 1);
            end
            if ~isempty(p.Results.y0)
                obj.solver_initialization.y0 = p.Results.y0;
            end

            [solution,stats] = obj.solve(obj.solver_initialization);
            varargout{1} = solution;
            obj.stats = stats;
        end

        function obj = parenAssign(obj,index_op,varargin)
            % nosnoc.error('invalid', 'Invalid operation');
            % TODO: there is no nosnoc in mpecopt - adjust errors messages
            error('mpecopt: Invalid operation.')
        end
        
        function obj = parenDelete(obj,index_op)
            % nosnoc.error('invalid', 'Invalid operation')
            % TODO: there is no nosnoc in mpecopt - adjust errors messages
            error('mpecopt: Invalid operation.')
        end

        function n = parenListLength(obj,index_op,ctx)
            n = 1;
        end
    end

    methods(Access=private)
        function process_solver_initialization(obj, solver_initialization)
            import casadi.*
            opts = obj.opts;
            mpec_casadi = obj.mpec_casadi;
            lpec_casadi = obj.lpec_casadi;
            dims = obj.dims;

            % check does a parameter exist:
            if ~isfield(solver_initialization,"p0")
                solver_initialization.p0 = [];
            end
            
            G_fun = mpec_casadi.G_fun;
            H_fun = mpec_casadi.H_fun;
            if opts.initial_comp_all_zero
                G_eval = zeros(dims.n_comp,1);
                H_eval = zeros(dims.n_comp,1);
            else
                G_eval = full(mpec_casadi.G_fun(solver_initialization.x0, solver_initialization.p0));
                H_eval = full(mpec_casadi.H_fun(solver_initialization.x0, solver_initialization.p0));
            end

            if opts.lift_complementarities_full
                solver_initialization.lbx = [solver_initialization.lbx;0*ones(dims.n_comp,1)];
                solver_initialization.ubx = [solver_initialization.ubx;inf*ones(dims.n_comp,1)];
                solver_initialization.lbg = [solver_initialization.lbg;0*ones(dims.n_comp,1)];
                solver_initialization.ubg = [solver_initialization.ubg;0*ones(dims.n_comp,1)];
                solver_initialization.x0 = [solver_initialization.x0;G_eval];
            else
                solver_initialization.lbx = [solver_initialization.lbx;0*ones(dims.n_lift_x1,1)];
                solver_initialization.ubx = [solver_initialization.ubx;inf*ones(dims.n_lift_x1 ,1)];
                solver_initialization.lbg = [solver_initialization.lbg;0*ones(dims.n_lift_x1 ,1)];
                solver_initialization.ubg = [solver_initialization.ubg;0*ones(dims.n_lift_x1 ,1)];
                solver_initialization.x0 = [solver_initialization.x0;G_eval(dims.ind_nonscalar_x1)];
            end

            if opts.lift_complementarities_full
                solver_initialization.lbx = [solver_initialization.lbx;0*ones(dims.n_comp,1)];
                solver_initialization.ubx = [solver_initialization.ubx;inf*ones(dims.n_comp,1)];
                solver_initialization.lbg = [solver_initialization.lbg;0*ones(dims.n_comp,1)];
                solver_initialization.ubg = [solver_initialization.ubg;0*ones(dims.n_comp,1)];
                solver_initialization.x0 = [solver_initialization.x0;H_eval];
            else
                solver_initialization.lbx = [solver_initialization.lbx;0*ones(dims.n_lift_x2,1)];
                solver_initialization.ubx = [solver_initialization.ubx;inf*ones(dims.n_lift_x2 ,1)];
                solver_initialization.lbg = [solver_initialization.lbg;0*ones(dims.n_lift_x2 ,1)];
                solver_initialization.ubg = [solver_initialization.ubg;0*ones(dims.n_lift_x2 ,1)];
                solver_initialization.x0 = [solver_initialization.x0;H_eval(dims.ind_nonscalar_x2)];
            end

            %% Split into equalites and inequalities
            % TODO@Anton?: Get rid of this unfold?
            x = mpec_casadi.x;
            g = mpec_casadi.g;
            x0 = mpec_casadi.x0;
            x1 = mpec_casadi.x1;
            x2 = mpec_casadi.x2;
            p = mpec_casadi.p;
            
            ind_g_eq = find(solver_initialization.lbg == solver_initialization.ubg);
            ind_g_ineq = find(solver_initialization.lbg < solver_initialization.ubg);

            ind_g_ineq_lb = find(solver_initialization.lbg >- inf & solver_initialization.lbg < solver_initialization.ubg);
            ind_g_ineq_ub = find(solver_initialization.ubg < inf & solver_initialization.lbg < solver_initialization.ubg);

            ind_x_lb = find(solver_initialization.lbx > -inf);
            ind_x_ub =  find(solver_initialization.ubx < inf);

            n_eq = length(ind_g_eq);
            n_g_ineq_ub = length(ind_g_ineq_ub);
            n_g_ineq_lb = length(ind_g_ineq_lb);

            n_ubx = length(ind_x_ub);
            n_lbx = length(ind_x_lb);

            lbx_reduced = solver_initialization.lbx(ind_x_lb);
            ubx_reduced = solver_initialization.ubx(ind_x_ub);
            % Generate casadi functions for objective and constraint function evaluations
            % Zero order
            g_eq = g(ind_g_eq)-solver_initialization.lbg(ind_g_eq);                    % g_eq = g - g_lb = 0
            g_ineq_ub = solver_initialization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
            g_ineq_lb = g(ind_g_ineq_lb)-solver_initialization.lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0
            g_ineq = [solver_initialization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
                g(ind_g_ineq_lb)-solver_initialization.lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
            n_ineq = size(g_ineq,1);
            % first-order constraint Jacobians
            if n_eq > 0
                nabla_g_eq = g_eq.jacobian(x);
            else
                nabla_g = [];
                nabla_g_eq = [];
            end
            if n_ineq > 0
                nabla_g_ineq_ub = g_ineq_ub.jacobian(x);
                nabla_g_ineq_lb = g_ineq_lb.jacobian(x);
            else
                nabla_g_ineq_ub = [];
                nabla_g_ineq_lb = [];
            end
            nabla_g_ineq = [nabla_g_ineq_ub; nabla_g_ineq_lb];

            %% Infeasiblity meausre in inf norm
            % all standard constraints
            h_eq = max(abs(g_eq));
            h_ineq_ub = max(min(g_ineq_ub,0));
            h_ineq_lb = max(min(g_ineq_lb,0));
            h_ubx = max(min(solver_initialization.ubx-x,0));
            h_lbx = max(min(x-solver_initialization.lbx,0));
            % Summary
            h_std = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);
            if dims.n_comp > 0
                if opts.comp_res_bilinear
                    h_comp= max(abs(x1).*abs(x2)); % here kappa =0.1 to have reduce of the value by factor of 10
                                                   % h_comp= sqrt(max(abs(min((x1),(x2))))); % here kappa =0.01 to reduce the value above by factor of 10
                else
                    h_comp= max(min(abs(x1),abs(x2))); % 
                    % h_comp = max(abs(min(x1,x2)));
                end
            else
                h_comp  = 0;
            end

            %% CasADi Functions for constraints infeasiblity
            mpec_casadi.h_std_fun  = Function('h_std_fun',{x,p},{h_std});
            mpec_casadi.h_comp_fun  = Function('h_comp_fun',{x,p},{h_comp});
            mpec_casadi.h_total_fun = Function('h_comp_fun',{x,p},{max(h_comp,h_std)});
            
            %% generate required functions
            % Zero order
            mpec_casadi.g_eq_fun =  Function('g_eq_fun',{x,p},{g_eq});
            mpec_casadi.g_ineq_ub_fun = Function('g_ineq_ub_fun',{x,p},{g_ineq_ub});
            mpec_casadi.g_ineq_lb_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq_lb});
            mpec_casadi.g_ineq_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq});
            % First order
            mpec_casadi.nabla_g_eq_fun = Function('nabla_g_eq_fun',{x,p},{nabla_g_eq});
            mpec_casadi.nabla_g_ineq_ub_fun  = Function('nabla_g_ineq_ub_fun',{x,p},{nabla_g_ineq_ub});
            mpec_casadi.nabla_g_ineq_lb_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq_lb});
            mpec_casadi.nabla_g_ineq_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq});

            %% lpec update
            % TODO(@anton) remove duplication
            sense_B = repmat('>',1,2*dims.n_comp);
            sense = [repmat('=',1,n_eq), repmat('>',1,n_ineq), sense_B];
            lpec_casadi.sense = sense;
            lpec_casadi.lbx = solver_initialization.lbx;
            lpec_casadi.ubx = solver_initialization.ubx;
            lpec_casadi.g_eq_fun = mpec_casadi.g_eq_fun;
            lpec_casadi.g_ineq_fun = mpec_casadi.g_ineq_fun;
            lpec_casadi.nabla_g_eq_fun = mpec_casadi.nabla_g_eq_fun;
            lpec_casadi.nabla_g_ineq_fun = mpec_casadi.nabla_g_ineq_fun;
            
            %% Update dims
            dims.n_eq = n_eq;
            dims.n_ineq = n_ineq;
            % index sets of general constraints
            dims.ind_g_eq = ind_g_eq;
            dims.ind_g_ineq = ind_g_ineq;
            dims.ind_g_ineq_lb = ind_g_ineq_lb;
            dims.ind_g_ineq_ub = ind_g_ineq_lb; % if the last two have same indicies then it is a two sided ineq;
            
            %% update structs
            obj.solver_initialization = solver_initialization;
            obj.mpec_casadi = mpec_casadi;
            obj.lpec_casadi = lpec_casadi;
            obj.dims = dims;
        end

        function create_mpec_functions(obj)
            import casadi.*
            mpec = obj.mpec;
            opts = obj.opts;

            % TODO(@anton) we need to choose one and stick to it.
            % TODO(@anton) Here we need to do interleaving again or do it in vdx.
            if isa(mpec, 'vdx.problems.Mpcc')
                x = mpec.w.sym;
                f = mpec.f;
                g = mpec.g.sym;
                G = mpec.G.sym;
                H = mpec.H.sym;
                p = mpec.p.sym;               
            else
                try
                    x = mpec.x;
                catch
                    x = mpec.w;
                end
                f = mpec.f;
                g = SX(mpec.g);
                G = mpec.G;
                H = mpec.H;

                if isfield(mpec,'p')
                    p = mpec.p;
                else
                    p = [];
                end
            end

            %% Edit complementarity constraints
            n_primal_non_lifted = length(x);
            n_comp = size(G,1);
            G_fun = Function('G_fun',{x,p},{G});
            H_fun = Function('H_fun',{x,p},{H});

            G_copy = G_fun(x,p);
            H_copy = H_fun(x,p);

            % check if the comps are expressions or just subvectors of w.
            if opts.lift_complementarities_full
                % define lift vairables
                x1 = SX.sym('x1',n_comp);
                % update x and init guess
                x = [x;x1];
                % lift
                g = [g;x1-G];
            else
                % lifting with only those that are not scaler
                % define lift vairables
                [ind_scalar,ind_nonscalar_x1, ind_map] = find_nonscalar(G,x);
                n_lift_x1 = length(ind_nonscalar_x1);
                if n_lift_x1 == 0
                    % TODO(@anton) Figure out what this does.
                    try
                        x.jacobian(G_copy);
                    catch
                        n_lift_x1 = length(G_copy);
                        ind_nonscalar_x1 = 1:n_lift_x1;
                        ind_scalar = [];
                    end
                end
                if n_lift_x1 > 0
                    x1_lift = SX.sym('x1_lift',n_lift_x1);
                    % x1 = [x(ind_scalar);x1_lift];
                    % update x and init guess
                    x = [x;x1_lift];
                    % lift
                    g = [g;x1_lift-G(ind_nonscalar_x1)];

                    x1 = G_copy;
                    x1(ind_nonscalar_x1) = x1_lift;
                else
                    x1 = G;
                end
            end

            if opts.lift_complementarities_full
                % define lift vairables
                x2 = SX.sym('x2',n_comp);
                % update x and init guess
                x = [x;x2];
                % lift
                g = [g;x2-H];
            else
                % lifting with only those that are not scaler
                [ind_scalar,ind_nonscalar_x2, ind_map] = find_nonscalar(H,x);
                n_lift_x2 = length(ind_nonscalar_x2);
                if n_lift_x2 == 0
                    % TODO(@anton) Figure out what this does.
                    try
                        x.jacobian(H_copy);
                    catch
                        n_lift_x2 = length(H_copy);
                        ind_nonscalar_x2 = 1:n_lift_x2;
                        ind_scalar = [];
                    end
                end
                if n_lift_x2 > 0
                    x2_lift = SX.sym('x2_lift',n_lift_x2);
                    % x2 = [x(ind_scalar);x2_lift];
                    % update x and init guess
                    x = [x;x2_lift];
                    % lift
                    g = [g;x2_lift-H(ind_nonscalar_x2)];

                    x2 = H_copy;
                    x2(ind_nonscalar_x2) = x2_lift;
                else
                    x2 = H;
                end
            end

            % find index set
            if n_comp > 0
                % TODO(@anton) do we actually need this here or can we calculate these "analytically"
                ind_x1 = [];
                ind_x2 = [];
                ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
                [ind_x1,~] = find(sparse(ind_x1_fun(zeros(size(x))))==1);
                ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
                [ind_x2,~] = find(sparse(ind_x2_fun(zeros(size(x))))==1);
                opts.nlp_is_mpec = 1; % TODO: check is this still used? (its purpose: if not an mpec, just make single nlp call without mpec machinery);
            else
                opts.nlp_is_mpec = 0;
                ind_x1 = [];
                ind_x2 = [];
            end

            n_primal = length(x);
            n_primal_x0 = n_primal - 2*n_comp; % primal variables excluding the complementarity variables;
            ind_x0 = [1:n_primal]';
            if opts.nlp_is_mpec
                ind_x0([ind_x1,ind_x2]) = [];
            end
            x0 = x(ind_x0); % Variables not involved in complementarity constraints.
            
            nabla_f = f.jacobian(x)';
            nabla_g = g.jacobian(x);
            %% CasADi functions for constraint evaluations and their derivaties
            % Zero order (objective and constraint function evaluations)
            mpec_casadi = struct;

            mpec_casadi.f = f;
            mpec_casadi.g = g;
            mpec_casadi.x = x;
            mpec_casadi.x1 = x1;
            mpec_casadi.x2 = x2;
            mpec_casadi.x0 = x0;
            mpec_casadi.p = p;

            mpec_casadi.f_fun = Function('f_fun',{x,p},{f});
            mpec_casadi.G_fun = G_fun;
            mpec_casadi.H_fun = H_fun;
            % First order (Gradients and Jacobian)
            mpec_casadi.nabla_f_fun = Function('nabla_f_fun',{x,p},{nabla_f});
            mpec_casadi.nabla_g_fun = Function('nabla_g_fun',{x,p},{nabla_g});

            %% Store some dimensions
            dims.n_slacks = 0; % in generla no slacks, except in feasiblity problems
            dims.ind_x0 = ind_x0;
            dims.ind_x1 = ind_x1;
            dims.ind_x2 = ind_x2;
            dims.n_primal = n_primal;
            dims.n_comp = n_comp;
            dims.n_primal_non_lifted = n_primal_non_lifted;
            dims.n_lift_x1 = n_lift_x1;
            dims.n_lift_x2 = n_lift_x2;
            dims.n_auxiliary = dims.n_comp; % number of binary variables in LPEC

            % indices for lifting
            dims.ind_nonscalar_x1 = ind_nonscalar_x1;
            dims.ind_nonscalar_x2 = ind_nonscalar_x2;

            %%  Main piece NLP solver (BNLP or TNLP)
            % TODO(@anton) actually build the scholtes relaxation here and store the scholtes indices.
            t_generate_nlp_solvers = tic;
            nlp = struct('x', x,'f', f,'g', g,'p',p);
            solver = nlpsol('solver', 'ipopt', nlp, opts.settings_casadi_nlp);
            stats.cpu_time_generate_nlp_solvers = toc(t_generate_nlp_solvers);
            mpec_casadi.solver = solver;

            %% Store structs
            obj.mpec_casadi = mpec_casadi;
            obj.dims = dims;
            obj.stats = stats;
        end

        function create_lpec_functions(obj)
            import casadi.*
            mpec_casadi = obj.mpec_casadi;
            dims = obj.dims;

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
            vtype = [repmat('C',1,dims.n_primal), repmat('B',1,dims.n_comp)];
            vtype_num = [zeros(1,dims.n_primal), ones(1,dims.n_comp)];

            lb_binary = 0*ones(dims.n_comp,1);
            ub_binary = 1*ones(dims.n_comp,1);

            lpec_casadi.A_lpec_fun = A_lpec_fun;
            lpec_casadi.b_res_fun = b_res_fun;
            lpec_casadi.vtype = vtype;
            lpec_casadi.vtype_num = vtype_num;
            lpec_casadi.lb_binary = lb_binary;
            lpec_casadi.ub_binary = ub_binary;
            % Copy mpec casadi functions needed for the linearization
            lpec_casadi.f_fun = mpec_casadi.f_fun;
            lpec_casadi.nabla_f_fun = mpec_casadi.nabla_f_fun;

            obj.lpec_casadi = lpec_casadi;
        end
    end
end
