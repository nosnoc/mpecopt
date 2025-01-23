% Copyright (c) 2024, Armin Nurkanović
% TODO: allow passing no lbg ubg and assume g empty in this case; and vice versa
classdef MPECOptimizerOptions < handle
    properties
        % General
        solver_name {mustBeTextScalar} = 'MPECopt';
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX';
        comp_res_bilinear (1,1) logical = true;  % If true: comp_res = max(x1.*x2), if false: comp_res = max(min(x1,x2)); the bilinear shrinks faster, e.g. x1 = 1e-3,x2 = 1e-3, bilinear = 1e-6, std 1e-3;
        lift_complementarities_full (1,1) logical = false; % If true brute force lifting x1-G = 0, if false, detect which components of G are already purely scalar. Affine expressions get lifted.

        %  --- Stopping criteria/tolernaces ---
        max_iter (1,1) double {mustBeInteger, mustBePositive} = 25;
        max_inner_iter (1,1) double {mustBeInteger, mustBePositive} = 6;
        
        tol (1,1) double {mustBeReal, mustBeNonnegative} = 1e-8;
        tol_B_stationarity (1,1) double {mustBeReal, mustBeNonnegative} = 1e-8;
        tol_feasibility (1,1) double {mustBeReal, mustBeNonnegative} = 1e-9;
        tol_B_stationarity_early_term (1,1) double {mustBeReal, mustBeNonnegative} = 1e-7;
        tol_active (1,1) double {mustBeReal, mustBeNonnegative} = 1e-9; % below this treshold a constaint is considered to be active
        
        rescale_large_objective_gradients (1,1) logical = true;
        allow_early_termination (1,1) logical = true; % if max time or max iters reached, accept problems with lower accuracy as sucess
        fast_B_stationarity_check (1,1) logical = false; % not supported anymore, esentially if x0 feasible, first solve a LPEC with an extreamly small TR radius 
        stop_if_S_stationary (1,1) logical = false; % S stationarity is sufficent for B

        % --- BNLP/TNLP settings----
        piece_nlp_strategy PieceNLPStartegy = "BNLP_integer";

        % Multiplier-based stationarity
        compute_tnlp_stationary_point (1,1) logical = true; % todo: reduce this and the one below settings to one
        resolve_tnlp_for_stationarity (1,1) logical = false;  % if the method is solving BNLPs, then resolve a TNLP to get correct MPEC multipliers
        plot_mpec_multipliers (1,1) logical = false;

        % ----- Phase I Settings -----
        initial_comp_all_zero(1,1) logical = false; % x1^0 = 0, x2^0 =0
        stop_if_nlp_infeasible (1,1) logical = true; % if an relaxed NLP is infeasible, declare the MPEC locally infeasible
        initalization_strategy (1,1) InitalizationStrategy = InitalizationStrategy.RelaxAndProject;
        max_recovery_iters (1,1) double {mustBeReal, mustBeInteger} = 12; % If the current BNLP is infeasible, try to solve a tighet relaxation of the MPEC for a better feasible BNLP guess;
        bnlp_projection_strategy (1,1) BNLPProjectionStrategy = BNLPProjectionStrategy.LPEC; % chose how to project x_k(tau) onto complementarity set
        project_guess_to_bounds (1,1) logical = false; % project point to comps and simple bounds (might be still infeasible for general constraints)
        

        relax_and_project_iters (1,1) double {mustBeReal, mustBeInteger} = 1;
        lpec_recovery_start (1,1) double {mustBeReal, mustBeNonnegative} = 0.2; % after how many % of failed recoverz iterations the LPEC recovery starts (meanse, TR of LPEC gets larger)
        lpec_solve_if_comp_feasible (1,1) logical = true; % if true dont try to solve lpecs if they are likely to be infeasible;
        relax_and_project_tighter_TR (1,1) logical = false; % if true take the maximum as max(max(x1,x2)), otherwise take it as rho_TR_phase_i
        relax_and_project_sigma0 (1,1) double {mustBeReal, mustBeNonnegative} = 1;
        relax_and_project_kappa (1,1) double {mustBeReal, mustBeNonnegative} = 0.1;
        relax_and_project_comp_tol (1,1) double {mustBeReal, mustBeNonnegative} = 1e-10;
        relax_and_project_comps_aggregated (1,1) logical = false;
        relax_and_project_consider_all_comps_in_lpec (1,1) logical = false;
        relax_and_project_homotopy_parameter_steering HomotopySteering = 'Direct';
        realx_and_project_scale_factor_rho_tr (1,1) double {mustBeReal, mustBeNonnegative} = 2; % in phase i rho_tr = scale_factor*rho_min_theortically
        fallback_to_intial_guess (1,1) logical = false;
        count_first_lpec_into_phase_i (1,1) logical = true; % if true, if first lpec in phase i verifies b-stationarity, delcare problem solved in phase i.

        % ------ Phase I feasiblity settings ----------
        feasibility_project_to_bounds  (1,1) logical = true;
        feasibility_s_lower_bound (1,1) double {mustBeReal} = 0.0; % lower bound on slack in feasbility I
        % trust_region_on_slacks  (1,1) logical = true; % Do the slack variables for the feasbility trasformation have a TR constraint?
        
        % ----- Phase II Settings -----
        rho_TR_min (1,1) double {mustBeReal, mustBePositive} = 1e-6; % TODO:should be slightly larger than the B_stat_tolerances, and should be in line with max inter iter
        rho_TR_max (1,1) double {mustBeReal, mustBePositive} = 1e6;
        rho_TR_phase_i_init (1,1) double {mustBeReal, mustBePositive} = 1e-1; % default starting value
        rho_TR_phase_ii_init (1,1) double {mustBeReal, mustBePositive} = 1e-3; % default starting value
        TR_increasing_factor(1,1) double {mustBeReal, mustBePositive} = 10;  % how much to incerase if sucess , 1.5
        TR_reducing_factor(1,1) double {mustBeReal, mustBePositive} = 0.1; % how much to reduce if fail, 0.75
        accept_last_inner_iter(1,1) logical = true; % if max number of inner iterations reached, accept step even if not sufficient decrees      
        smaller_tr_in_phase_ii (1,1) logical = false; % TODO: align with delta tr abooveonce feasiblity is reached, solve lpecs with a very small tr.
        reset_TR_radius (1,1) logical = true;
    
        % ----- LPEC solver settings ----
        consider_all_complementarities_in_lpec(1,1) logical = true; % if true, take all comps into the lpec, if false take only those active at x^k
        reduced_lpec_via_fixed_integers(1,1) logical = true; % if true, do not reduce the matrix B_lpec, but just fix the known binary variables; - FALSE BROKEN
        warm_start_lpec_solver(1,1) logical = true;
        tighten_bounds_in_lpec(1,1) logical = false; %TODO: remove
        BigM(1,1) double {mustBeReal, mustBePositive} = 1e2;
        constant_lpec_TR_radius(1,1) logical = false;  % TODO: remove more or less useles option
        rho_TR_lpec(1,1) double {mustBeReal, mustBePositive} = 1e1; % the TR radius in the lpec if constant.

        %  LPEC fallback strategy
        use_lpec_fallback_strategy (1,1) logical = false; 
        use_lpec_fallback_strategy_phase_i (1,1) logical = false; % if infeasible make TR larger, if node limit reached, make TR radius smaller
        use_lpec_fallback_strategy_phase_ii (1,1) logical = false; % if infeasible make TR larger, if node limit reached, make TR radius smaller
        take_TR_from_fallback_strategy (1,1) logical = false; % take for trust region radius the one from the last triggered fall back strategy (usually a large number)
        fallback_TR_increasing_factor (1,1) double {mustBeReal, mustBePositive} = 10;  %
        fallback_TR_reducing_factor (1,1) double {mustBeReal, mustBePositive} = 0.1; %
        N_attemps_infeasible (1,1) double {mustBeInteger, mustBeNonnegative} =  5;
        N_attemps_node_limit (1,1) double {mustBeInteger, mustBeNonnegative} = 2;

        % ---------NLP solver settings------
        compute_bnlp_step (1,1) logical = true;
        use_initial_nlp_constraint_ordering(1,1) logical = true; %  true =  use in nlp  lbg <= g <ubg, false use reordered equality/inequality like in LPEC
        initalize_bnlp_with_lpec_solution(1,1) logical = false;
        nlp_is_mpec(1,1) logical = true;

        % ----- subsolver settings ----------
        settings_casadi_nlp % NLP solver settings
        settings_lpec %LPEC solver settings

        % ---- Debug ---------
        plot_lpec_iterate (1,1) logical = false;
        debug_mode_on (1,1) logical = false;
        stop_if_step_rejected (1,1) logical = false;  
        reset_lpec_objective (1,1) logical = true;

        % ---- Print -----
        verbose_solver(1,1) logical = true;
        verbose_summary(1,1) logical = true;
        verbose_extended_summary(1,1) logical = false;
    end

    methods
        function obj = MPECOptimizerOptions()
            obj.settings_lpec = LPECSolverOptions();
            % obj.settings_lpec.solver_name = 'lpec_solver';
            % obj.settings_lpec.lpec_solver= LpecSolver.Gurobi;
            % obj.settings_lpec.max_nodes = 2e3;
            obj.settings_lpec.stop_lpec_at_feasible = false;
            obj.settings_lpec.trust_region_on_slacks = false; % trust region for slack variables in fesability problem;
            % obj.settings_lpec.rel_tol= 1e-6;
            % obj.settings_lpec.abs_tol= 1e-8;
        
            default_tol = 1e-12;
            obj.settings_casadi_nlp.ipopt.print_level = 0;
            obj.settings_casadi_nlp.print_time = 0;
            obj.settings_casadi_nlp.ipopt.sb = 'yes';
            obj.settings_casadi_nlp.verbose = false;
            obj.settings_casadi_nlp.ipopt.max_iter = 3000;
            obj.settings_casadi_nlp.ipopt.bound_relax_factor = 0;
            % obj.settings_casadi_nlp.ipopt.bound_relax_factor = 1e-16;
            % obj.settings_casadi_nlp.ipopt.honor_original_bounds = 'yes';
            obj.settings_casadi_nlp.ipopt.tol = default_tol;
            obj.settings_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.settings_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.settings_casadi_nlp.ipopt.compl_inf_tol = default_tol;
            obj.settings_casadi_nlp.ipopt.acceptable_tol = 1e-9;
            obj.settings_casadi_nlp.ipopt.mu_strategy = 'adaptive';
            obj.settings_casadi_nlp.ipopt.mu_oracle = 'quality-function';
            obj.settings_casadi_nlp.ipopt.warm_start_init_point = 'yes';
            % obj.settings_casadi_nlp.ipopt.warm_start_bound_push = 1e-15;
            % obj.settings_casadi_nlp.ipopt.warm_start_bound_frac = 1e-15;
            % obj.settings_casadi_nlp.ipopt.warm_start_bound_push = 1e-6;
            % obj.settings_casadi_nlp.ipopt.warm_start_bound_frac = 1e-6;
            obj.settings_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
            obj.settings_casadi_nlp.ipopt.linear_solver = 'mumps'; % mumps, ma27, ma57
            obj.settings_casadi_nlp.ipopt.fixed_variable_treatment  = 'make_parameter';  % make_parameter  make_constraint relax_bounds
            
            % obj.opts_casadi_nlp.snopt = struct();
            % obj.opts_casadi_nlp.worhp = struct();
            % obj.opts_casadi_nlp.uno = struct();
            % obj.p_val = [obj.sigma_0];
        end

        function [] = preprocess(obj)
            %             import casadi.*
        end

    end
end
