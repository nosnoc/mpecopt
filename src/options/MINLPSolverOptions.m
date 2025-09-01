% Copyright (c) 2025, Armin Nurkanović
classdef MINLPSolverOptions< handle
    properties
        % General
        solver_name {mustBeTextScalar} = 'mpec_minlp';
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX';
        comp_res_bilinear(1,1) logical = true;  % if true comp_res = max(x1.*x2), if false, comp_res = max(min(x1,x2)); the bliniear shrinks faster, e.g. x1 = 1e-3,x2 = 1e-3, bilinear = 1e-6, std 1e-3;

        % Stopping criteria/tolernaces
        max_iter(1,1) double {mustBeInteger, mustBePositive} = 15;
        tol(1,1) double {mustBeReal, mustBeNonnegative} = 1e-9;
        tol_active(1,1) double {mustBeReal, mustBeNonnegative} = 1e-10; % below this treshold a constraint is considered to be active
        comp_tol(1,1) double {mustBeReal, mustBeNonnegative} = 1e-9;
        plot_mpec_multipliers(1,1) logical = false;
        initial_comp_all_zero(1,1) logical = false;

        success_only_if_s_stationary(1,1) logical = false;
        compute_tnlp_stationary_point(1,1) logical = true;
        check_B_stationarity(1,1) logical = true;
        tol_B_stationarity(1,1) double {mustBeReal, mustBeNonnegative} = 1e-8;
        consider_all_complementarities_in_lpec(1,1) logical = true; % for B stat check
        rescale_large_objective_gradients (1,1) logical = true;

        lift_complementarities(1,1) logical = true; % smart lifting, dont introduce slack variable for those that are already scalar;
        lift_complementarities_full(1,1) logical = false; % brute force lifting without checking are some parts of G or H already scalar;
        aggregate_comps(1,1) logical = false;

        % LPEC settings (needed in B stationarity detection)
        problem_is_lpec(1,1) logical = false;
        update_comp_lower_bounds(1,1) logical = true; % if true set lbx and ubx of comp constraints to zero
        BigM(1,1) double {mustBeReal, mustBeNonnegative} = 1e4; % LPEC
        BigM_minlp(1,1) double {mustBeReal, mustBeNonnegative} = 1e4; % MINLP

        % verbose
        verbose_solver(1,1) logical = true;
        verbose_summary(1,1) logical = true;

        % NLP solver settings
        settings_casadi_nlp

    end

    methods
        function obj = MINLPSolverOptions()
            default_tol = 1e-10;
            % --------- MINLP Specific settings
            % bonmin: https://www.coin-or.org/Bonmin/option_pages/options_list_bonmin.html
            obj.settings_casadi_nlp.error_on_fail = false;
            obj.settings_casadi_nlp.bonmin.algorithm = 'B-BB';
            obj.settings_casadi_nlp.bonmin.warm_start = 'interior_point';
            obj.settings_casadi_nlp.bonmin.node_limit = 75;
            obj.settings_casadi_nlp.bonmin.time_limit = 600;

            % options for the algorithm
            % B-BB % simple branch-and-bound algorithm,
            % B-OA % OA Decomposition algorithm,
            % B-QG % Quesada and Grossmann branch-and-cut algorithm,
            % B-Hyb % hybrid outer approximation based branch-and-cut,
            % B-Ecp % ECP cuts based branch-and-cut a la FilMINT.
            % B-iFP % Iterated Feasibility Pump for MINLP.

            % ------ NLP Solver settings for ipopt ---------------
            % ipopt: https://coin-or.github.io/Ipopt/OPTIONS.html
            obj.settings_casadi_nlp.bonmin.print_level = 0;
            obj.settings_casadi_nlp.bonmin.print_level = 0;
            obj.settings_casadi_nlp.print_time = 0;
            obj.settings_casadi_nlp.verbose = false;
            obj.settings_casadi_nlp.bonmin.max_iter = 3000;
            % obj.settings_casadi_nlp.bonmin.bound_relax_factor = 0;
            obj.settings_casadi_nlp.bonmin.bound_relax_factor = default_tol;
            obj.settings_casadi_nlp.bonmin.honor_original_bounds = 'yes';
            obj.settings_casadi_nlp.bonmin.tol = default_tol;
            obj.settings_casadi_nlp.bonmin.dual_inf_tol = default_tol;
            obj.settings_casadi_nlp.bonmin.dual_inf_tol = default_tol;
            obj.settings_casadi_nlp.bonmin.compl_inf_tol = default_tol;
            obj.settings_casadi_nlp.bonmin.acceptable_tol = 1e-9;
            obj.settings_casadi_nlp.bonmin.mu_strategy = 'adaptive';
            obj.settings_casadi_nlp.bonmin.mu_oracle = 'quality-function';
            % obj.settings_casadi_nlp.bonmin.warm_start_init_point = 'yes';
            % obj.settings_casadi_nlp.bonmin.warm_start_entire_iterate = 'yes';
            obj.settings_casadi_nlp.bonmin.linear_solver = 'ma27'; % 'mumps'; ma57

            obj.settings_casadi_nlp.detect_simple_bounds = true;
            obj.settings_casadi_nlp.bonmin.fixed_variable_treatment  = 'relax_bounds';  % make_parameter  make_constraint relax_bounds

            % Some BonminSettings
            % obj.opts_casadi_nlp.snopt = struct();
            % obj.opts_casadi_nlp.worhp = struct();
            % obj.opts_casadi_nlp.uno = struct();
            % obj.p_val = [obj.sigma_0];
        end

    end
end
