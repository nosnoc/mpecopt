% Copyright (c) 2024, Armin Nurkanović
classdef LPECSolverOptions< handle
    properties
        % General
        solver_name {mustBeTextScalar} = 'lpec_solver';

        % MILP Solver for LPEC settings
        lpec_solver(1,1) LpecSolver = LpecSolver.Gurobi
        max_nodes(1,1) double {mustBeInteger, mustBePositive} = 5e2;
        max_time(1,1) double {mustBeReal, mustBePositive} = 2e2; % ~3 min time out;
        cutoff(1,1) double {mustBeReal} = 10;
        rel_tol(1,1) double {mustBeReal, mustBeNonempty} = 1e-9;
        abs_tol(1,1) double {mustBeReal, mustBeNonempty} = 1e-9;
        solve_lpec_with_cutoff (1,1) logical = false;
        stop_lpec_at_feasible (1,1) logical = false;
        trust_region_on_slacks  (1,1) logical = false; % Do the slack variables for the feasbility trasformation have a TR constraint?
        homotopy_solver_settings 

        % settings for custom nlp based lpec methods (Scholtes, Ell1, Ell_inf)
        max_iter(1,1) double {mustBeInteger, mustBePositive} = 25;
        tol(1,1) double {mustBeReal, mustBeNonnegative} = 1e-9; 
        aggregate_comps(1,1) logical = false;
        sigma0(1,1) double {mustBeReal, mustBeNonnegative} = 1;
        kappa(1,1) double {mustBeReal, mustBePositive} = 0.1;

    end

    methods
        function obj = LPECSolverOptions()
            obj.homotopy_solver_settings = HomotopySolverOptions();
            obj.homotopy_solver_settings.verbose_solver = 1;
            obj.homotopy_solver_settings.verbose_summary = 0;
            obj.homotopy_solver_settings.compute_tnlp_stationary_point = 0;
            obj.homotopy_solver_settings.check_B_stationarity = 0;
            obj.homotopy_solver_settings.lift_complementarities = 0;
            obj.homotopy_solver_settings.lift_complementarities_full = 0;
            obj.homotopy_solver_settings.update_comp_lower_bounds = 1;
            obj.homotopy_solver_settings.problem_is_lpec = true;
            obj.homotopy_solver_settings.comp_tol = 1e-9;
            obj.homotopy_solver_settings.tol = 1e-9;
            obj.homotopy_solver_settings.max_iter = 11;
            obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.acceptable_tol = 1e-9;
            obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.max_iter = 1500;
            obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.max_wall_time = 120;
            %  obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.linear_solver = 'ma27';
            % obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.linear_solver= 'mumps';
            % obj.homotopy_solver_settings.settings_casadi_nlp.ipopt.linear_solver= 'mumps';
            
            obj.homotopy_solver_settings.solver_name = 'lpec_solver_scholtes_ext';
        end

        function [] = preprocess(obj)
%             import casadi.*
        end

    end
end
