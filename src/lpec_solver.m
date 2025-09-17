function  [results,stats] = lpec_solver(lpec,settings)
%LPEC_SOLVER Solves Linear Programs with Complementarity Constraints
%
%   [RESULTS, STATS] = LPEC_SOLVER(LPEC, SETTINGS) solves LPECs of the form:
%       min_d   f' * d
%       s.t.    A_eq * d + b_eq = 0
%               A_ineq * d + b_ineq >= 0  
%               lb <= d + x_lin <= ub
%               0 <= d_1 + x_lin_1 ⊥ d_2 + x_lin_2 >= 0
%               -rho_TR <= d <= rho_TR
%
%   METHODS: Big-M with MILP solvers (Gurobi, HiGHS, MATLAB) or continuous 
%   reformulations (Scholtes, L1/L∞ penalty) with Ipopt via homotopy.
%
%   INPUTS:
%   lpec - Structure with fields:
%       .x_lin, .d_lpec, .y_lpec  - Linearization point and initial guesses
%       .rho_TR, .lb, .ub         - Trust region radius and bounds  
%       .f, .A_eq, .b_eq          - Objective and equality constraints
%       .A_ineq, .b_ineq          - Inequality constraints
%       .A_lpec, .b_lpec          - Complementarity constraint data
%       .dims, .sense, .vtype     - Problem dimensions and variable types
%   
%   settings - Solver options structure
%
%   OUTPUTS:
%   results - Solution: .d_lpec (cont. variables), y_lpec (binaries encoding active set), .f_opt (lpec objective)
%   stats   - Statistics: .solve_time, .iterations, .nodecount, .itercount, .gap
%%
import casadi.*
%% init some stats:
stats.itercount = 0; % number simplex iters
stats.baritercount = 0; % nuber of barrier iters
stats.nodecount = 0; % number of nodes in BnB
stats.gap = inf; % integer gap
stats.optimal_solution_found = false; % B-stationarity can be checked only if the LPEC is solved to optimality, in other cases it is sufficent to have a feasible point.

%% Prepare LPEC
% add the boundso of the inaries
if isempty(lpec.lb_binary) || isempty(lpec.ub_binary)
    lb = [lpec.lb;0*ones(lpec.dims.n_auxiliary,1)];
    ub = [lpec.ub;1*ones(lpec.dims.n_auxiliary,1)];
else
    lb = [lpec.lb;lpec.lb_binary];
    ub = [lpec.ub;lpec.ub_binary];
end

%% x0 lbounds
if lpec.dims.n_slacks > 0
    % Dont put TR bounds on feasiblity slacks
    temp = lpec.dims.ind_x0;
    if ~settings.trust_region_on_slacks
        lpec.dims.ind_x0(1:lpec.dims.n_slacks) = [];
        lb(1:lpec.dims.n_slacks) = lb(1:lpec.dims.n_slacks)-lpec.x_lin(1:lpec.dims.n_slacks);
        lb(lpec.dims.ind_x0) = max([lb(lpec.dims.ind_x0)-lpec.x_lin(lpec.dims.ind_x0),-lpec.rho_TR*ones(size(lpec.dims.ind_x0))]');
        lpec.dims.ind_x0 = temp;
    end
else
    lb(lpec.dims.ind_x0) = max([lb(lpec.dims.ind_x0)-lpec.x_lin(lpec.dims.ind_x0),-lpec.rho_TR*ones(size(lpec.dims.ind_x0))]');
end


%% x1 and x2  lbounds (complementarity variables)
lb(lpec.dims.ind_x1) = max([lb(lpec.dims.ind_x1)-lpec.x_lin(lpec.dims.ind_x1),-lpec.rho_TR*ones(lpec.dims.n_comp,1),-lpec.x_lin(lpec.dims.ind_x1)],[],2);
lb(lpec.dims.ind_x2) = max([lb(lpec.dims.ind_x2)-lpec.x_lin(lpec.dims.ind_x2),-lpec.rho_TR*ones(lpec.dims.n_comp,1),-lpec.x_lin(lpec.dims.ind_x2)],[],2);

% Uper bound on all
if settings.trust_region_on_slacks
    ub(1:lpec.dims.n_primal) = min(ub(1:lpec.dims.n_primal)-lpec.x_lin(1:lpec.dims.n_primal),lpec.rho_TR);
else
    ub(1:lpec.dims.n_slacks) = ub(1:lpec.dims.n_slacks)-lpec.x_lin(1:lpec.dims.n_slacks);
    ub(lpec.dims.n_slacks+1:lpec.dims.n_primal) = min(ub(lpec.dims.n_slacks+1:lpec.dims.n_primal)-lpec.x_lin(lpec.dims.n_slacks+1:lpec.dims.n_primal),lpec.rho_TR);
end


%% intial guess
if lpec.dims.n_auxiliary < lpec.dims.n_comp
    y_lpec = round(rand(lpec.dims.n_auxiliary,1));
end

if any(isnan(lpec.d_lpec)) || any(isnan(lpec.y_lpec))
    x0 = [zeros(lpec.dims.n_primal,1); zeros(lpec.dims.n_auxiliary,1)];
else
    x0 = [lpec.d_lpec; lpec.y_lpec];
end

f_lpec = [lpec.f;zeros(lpec.dims.n_auxiliary,1)]; % augmented objective

%% Set up problem specific to solver
switch settings.lpec_solver
    case "Gurobi"
        A_gurobi = sparse([[lpec.A_eq, zeros(lpec.dims.n_eq, lpec.dims.n_auxiliary)];...
            [lpec.A_ineq,zeros(lpec.dims.n_ineq, lpec.dims.n_auxiliary)];...
            lpec.A_lpec]);
        b_gurobi = [-lpec.b_eq; -lpec.b_ineq; lpec.b_lpec];

        model.A = A_gurobi;
        model.rhs = b_gurobi;
        model.sense = lpec.sense;
        model.obj = f_lpec;
        model.vtype = lpec.vtype;
        model.modelsense = 'min';
        model.lb = lb;
        model.ub = ub;
        params.outputflag = 0;
        params.NodeLimit = settings.max_nodes;
        % accuracy
        params.FeasibilityTol = 1e-9;
        params.BarConvTol = 1e-9;
        params.IntFeasTol = 1e-9;
        params.TimeLimit = settings.max_time;
        params.OptimalityTol = 1e-9;
        if settings.stop_lpec_at_feasible && settings.is_in_phase_i
            % terminate phase i lpecs at a feasible point
            % params.MIPGap = 1;
            params.SolutionLimit = 1;
        end
        if settings.stop_lpec_at_descent && ~settings.is_in_phase_i
            % terminate phase ii lpecs at a dscent direction
            params.BestObjStop = -0.1*norm(f_lpec)*lpec.rho_TR;
        end
        % params.ObjScale = -0.5;     % https://www.gurobi.com/documentation/current/refman/objscale.html#parameter:ObjScale
        % params.ScaleFlag=0; % default -1, % https://www.gurobi.com/documentation/current/refman/scaleflag.html#parameter:ScaleFlag
        if settings.solve_lpec_with_cutoff
            params.cutoff = settings.cutoff;
        end
        model.start = x0;
    case {"Highs","Matlab"}
        % create MILP
        A_ineq_matlab = ([-[lpec.A_ineq,zeros(lpec.dims.n_ineq, lpec.dims.n_auxiliary)];...
            -lpec.A_lpec]);
        b_ineq_matlab = [lpec.b_ineq; -lpec.b_lpec];
        intcon = find(lpec.vtype_num);
        A_eq_matlab = [lpec.A_eq, zeros(lpec.dims.n_eq, lpec.dims.n_auxiliary)];
        b_eq_matlab = -lpec.b_eq;
        model.f = f_lpec;
        model.intcon = intcon;
        model.A = A_ineq_matlab;
        model.b = b_ineq_matlab;
        model.Aeq = A_eq_matlab;
        model.beq = b_eq_matlab;
        model.lb = lb;
        model.ub = ub;
        % model.X0 = x0;
        options = optimoptions('intlinprog','Display','off');
        % options = optimoptions('intlinprog');
        if settings.solve_lpec_with_cutoff
            options.ObjectiveCutOff  = settings.cutoff;
        end
        if settings.stop_lpec_at_feasible && settings.is_in_phase_i
            options.RelativeGapTolerance = 1; % stop at first integer feasible sol
        end

        % GapTolerance = 1e-16;
        % GapTolerance  = 1e-4;
        options.MaxNodes = settings.max_nodes;
        options.ConstraintTolerance = 1e-9;
        % options.AbsoluteGapTolerance = GapTolerance;
        % options.RelativeGapTolerance = GapTolerance;
        options.MaxTime = settings.max_time;
        model.options = options;
        model.solver = 'intlinprog';
        if isequal(settings.lpec_solver,"Highs")
            model.Algorithm = 'highs';
        else
            model.Algorithm = 'legacy';
        end
    case "Highs_casadi"
        % See for options: https://ergo-code.github.io/HiGHS/dev/options/definitions/
        A_highs = sparse([[lpec.A_eq, zeros(lpec.dims.n_eq, lpec.dims.n_auxiliary)];...
            [lpec.A_ineq,zeros(lpec.dims.n_ineq, lpec.dims.n_auxiliary)];...
            lpec.A_lpec]);
        lbA_highs = [-lpec.b_eq; -lpec.b_ineq; lpec.b_lpec];
        ubA_highs = [-lpec.b_eq; inf(lpec.dims.n_ineq+2*lpec.dims.n_auxiliary,1)];
        c_highs = f_lpec;
        discrete = lpec.vtype;

        lp.a = DM(A_highs).sparsity();
        highs_opts = struct;
        highs_opts.discrete = lpec.vtype_num;
        highs_opts.highs.log_to_console = false;
        % highs_opts.highs.simplex_strategy = 4;
        % highs_opts.highs.error_on_fail = false;
        if settings.stop_lpec_at_feasible && settings.is_in_phase_i
             highs_opts.highs.mip_rel_gap = 1;
        end

        % highs_opts.highs.mip_max_nodes = settings.max_nodes;
        % highs_opts.highs.time_limit = settings.max_time;
        % highs_opts.highs.optimality_tolerance = 1e-9;
        % highs_opts.highs.mip_feasibility_tolerance = 1e-9;
        % highs_opts.highs.kkt_tolerance = 1e-9;
        % highs_opts.highs.primal_feasibility_tolerance = 1e-9;
        % highs_opts.highs.dual_feasibility_tolerance = 1e-9;        
        highs_opts.highs.mip_heuristic_effort = 0.1;

        % highs_opts.highs.simplex_strategy = 4;
        % highs_opts.error_on_fail = false;
        lpsol = conic('lp', 'highs', lp, highs_opts);

    case "Projected_Gradient"
        % preparation - project into feasible set;
        x_lin_0 = lpec.x_lin(lpec.dims.ind_x0);%
        x_lin_1 = lpec.x_lin(lpec.dims.ind_x1);%
        x_lin_2 = lpec.x_lin(lpec.dims.ind_x2);%

        ind_x0_infeasible = x_lin_0<lpec.lb(lpec.dims.ind_x0) | x_lin_0>lpec.ub(lpec.dims.ind_x0);
        % if not in bounds put into middle, if bounds exist; if upper or lower bound -inf, set to zero;
        x_lin_0(ind_x0_infeasible) = 0.5*(max(lpec.lb(ind_x0_infeasible),-1e3)+min(lpec.ub(ind_x0_infeasible),1e3)); % mean between reduced bounds
        x_lin_1 = max(0,x_lin_1);
        x_lin_2 = max(0,x_lin_2);
        x_lin_1(x_lin_1<=x_lin_2) = 0;
        x_lin_2(x_lin_1>x_lin_2) = 0; % cf. Section 2.3 in Kirches
    case {"Reg","Ell_1", "Ell_inf", 'Nlp'}
        % Shared code for all regularization-based LPEC solevers;
        d_lpec_sym = SX.sym('d_lpec_sym',lpec.dims.n_primal);
        f_lpec = lpec.f'*d_lpec_sym;
        G_lpec = d_lpec_sym(lpec.dims.ind_x1)+lpec.x_lin(lpec.dims.ind_x1);
        H_lpec = d_lpec_sym(lpec.dims.ind_x2)+lpec.x_lin(lpec.dims.ind_x2);
        x0 = lpec.d_lpec;
        % create lpec
        g = [];lbg = []; ubg = [];
        if lpec.dims.n_eq > 0
            g = [g; lpec.A_eq*d_lpec_sym + lpec.b_eq];
            lbg = [lbg; zeros(lpec.dims.n_eq,1)];
            ubg = [ubg; zeros(lpec.dims.n_eq,1)];
        end
        if lpec.dims.n_ineq > 0
            g = [g; lpec.A_ineq*d_lpec_sym + lpec.b_ineq];
            lbg = [lbg; zeros(lpec.dims.n_ineq,1)];
            ubg = [ubg; inf*ones(lpec.dims.n_ineq,1)];
        end
        lb = [max([-lpec.rho_TR*ones(lpec.dims.n_primal-2*lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x0)-lpec.x_lin(lpec.dims.ind_x0)],[],2);...
            max([-lpec.rho_TR*ones(lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x1)-lpec.x_lin(lpec.dims.ind_x1),-lpec.x_lin(lpec.dims.ind_x1)],[],2);...
            max([-lpec.rho_TR*ones(lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x2)-lpec.x_lin(lpec.dims.ind_x2),-lpec.x_lin(lpec.dims.ind_x2)],[],2)];
        ub = min([lpec.rho_TR*ones(lpec.dims.n_primal), lpec.ub - lpec.x_lin],[],2);

        h_comp = max(min(abs(G_lpec),abs(H_lpec)));
        h_comp_fun = Function('h_comp_fun',{d_lpec_sym},{h_comp});
        sigma_lpec = SX.sym('sigma_lpec',1);
        sigma_k_lpec = settings.sigma0;
        lpec_homotopy = struct('x',d_lpec_sym,'f',f_lpec,'g', g,'G',G_lpec,'H',H_lpec);
        solver_initalization = struct('x0', x0,'lbx',lb,'ubx',ub,'lbg',lbg,'ubg',ubg);
    case "LCLP"
        L = zeros([length(lpec.dims.ind_x1), lpec.dims.n_primal]);
        R = zeros([length(lpec.dims.ind_x1), lpec.dims.n_primal]);
        for ii=1:length(lpec.dims.ind_x1)
            L(ii, lpec.dims.ind_x1(ii)) = 1;
            R(ii, lpec.dims.ind_x2(ii)) = 1;
        end
        lb = [max([-lpec.rho_TR*ones(lpec.dims.n_primal-2*lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x0)-lpec.x_lin(lpec.dims.ind_x0)],[],2);...
            max([-lpec.rho_TR*ones(lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x1)-lpec.x_lin(lpec.dims.ind_x1),-lpec.x_lin(lpec.dims.ind_x1)],[],2);...
            max([-lpec.rho_TR*ones(lpec.dims.n_comp,1), lpec.lb(lpec.dims.ind_x2)-lpec.x_lin(lpec.dims.ind_x2),-lpec.x_lin(lpec.dims.ind_x2)],[],2)];
        ub = min([lpec.rho_TR*ones(lpec.dims.n_primal), lpec.ub - lpec.x_lin],[],2);
end

%% Solve Problem
switch settings.lpec_solver
    case "Gurobi"
        try
            gurobi_time = tic;
            result_gurobi = gurobi(model, params);
            stats.itercount = result_gurobi.itercount;
            stats.baritercount = result_gurobi.baritercount;
            stats.nodecount = result_gurobi.nodecount;
            stats.gap = result_gurobi.mipgap;
            cpu_time_gurobi = toc(gurobi_time);
            % todo; add node and itter coumt, and add them tho phase i iter
            % and node, 
        catch
            model;
            % keyboard;
            % result_gurobi = [];
            result_gurobi.status = 'error';
            % results.d_lpec = lpec.d_lpec*nan;
            % results.y_lpec = lpec.y_lpec*nan;
            % results.f_opt = nan;
            % stats.lpec_solution_exists = false;
            result_gurobi.nodecount = nan;
            result_gurobi.runtime = nan;
            result_gurobi.mipgap = nan;

            cpu_time_gurobi = nan;
        end

        if isequal(result_gurobi.status,'OPTIMAL')
            stats.optimal_solution_found = true;
        end

        if (isequal(result_gurobi.status,'OPTIMAL') || isequal(result_gurobi.status,'NODE_LIMIT')|| isequal(result_gurobi.status,'USER_OBJ_LIMIT') || isequal(result_gurobi.status,'SOLUTION_LIMIT')) && isfield(result_gurobi,'x')
            results.d_lpec = result_gurobi.x(1:lpec.dims.n_primal);
            results.y_lpec = result_gurobi.x(end-lpec.dims.n_auxiliary+1:end);
            results.f_opt = result_gurobi.objval;
            stats.lpec_solution_exists = true;
            % stats.nodecount = result_gurobi.nodecount;
        else
            results.d_lpec = lpec.d_lpec*nan;
            results.y_lpec = lpec.y_lpec*nan;
            results.f_opt = nan;
            stats.lpec_solution_exists = false;
            % stats.nodecount = result_gurobi.nodecount;
        end
        stats.solver_message = result_gurobi.status;
        % stats.cpu_time = result_gurobi.runtime;
        stats.cpu_time  = cpu_time_gurobi;
    case {"Highs", "Matlab"}
        try
            intlinprog_time = tic;
            % [x,f_opt,statsu,output] = intlinprog(model);
            [x,f_opt,status,output] = intlinprog(f_lpec,intcon,A_ineq_matlab,b_ineq_matlab,A_eq_matlab,b_eq_matlab,lb,ub,x0,options);
            cpu_time = toc(intlinprog_time);
        catch
            model;
            % keyboard;
            result_gurobi = [];
        end
        switch status
            case {1,2,3}
                results.d_lpec = x(lpec.vtype_num==0);
                results.y_lpec = round(x(lpec.vtype_num==1)); % sometimes numerical errors in binary in highs;
                results.f_opt = f_opt;
                stats.lpec_solution_exists = true;
                stats.nodecount = output.numnodes;
                stats.gap = output.relativegap;
                stats.solver_message  = 'OPTIMAL';
                if status == 2
                    stats.solver_message  = 'NODE_LIMIT';                  % node limit but solution exists
                else
                    stats.solver_message  = 'OPTIMAL';
                end
                if status == 1
                    stats.optimal_solution_found = true;
                end
            case {0,-2,-3,-9}
                results.d_lpec = lpec.d_lpec*nan;
                results.y_lpec = lpec.y_lpec*nan;
                results.f_opt = nan;
                stats.lpec_solution_exists = false;
                stats.nodecount = output.numnodes;
                stats.gap = nan;
                if status == 0
                    stats.solver_message  = 'NODE_LIMIT';
                else
                    stats.solver_message  = 'INFEASIBLE';
                end

        end
        stats.solver_message_extended = output.message;
        stats.cpu_time = cpu_time;
    case "Highs_casadi"
        higs_casadi_time = tic;
        try
            r = lpsol('g', c_highs, 'a', A_highs, 'lbx', lb, 'ubx', ub, 'lba', lbA_highs, 'uba', ubA_highs);
            highs_success = strcmp(lpsol.stats.return_status, 'Optimal');
        catch
            highs_success  = false;
        end
        cpu_time = toc(higs_casadi_time);

        if highs_success
            x = full(r.x);
            results.d_lpec = x(lpec.vtype_num==0);
            results.y_lpec = round(x(lpec.vtype_num==1));
            results.f_opt = full(r.cost);
            stats.lpec_solution_exists = true;
            stats.success = highs_success;

            if strcmp(lpsol.stats.return_status, 'Optimal');
                stats.optimal_solution_found = true;
            end
        else
            results.d_lpec = lpec.d_lpec*nan;
            results.y_lpec = lpec.y_lpec*nan;
            results.f_opt = nan;
            stats.lpec_solution_exists = false;
        end
        stats.nodecount = lpsol.stats.n_call_solver;
        stats.itercount = lpsol.stats.simplex_iteration_count;
        stats.baritercount = lpsol.stats.ipm_iteration_count;
        stats.gap = lpsol.stats.mip_gap;
        % stats.cpu_time =  lpsol.stats.t_wall_solver;
        stats.cpu_time  = cpu_time;
        

        stats.solver_message =  lpsol.stats.unified_return_status;
    case {'Reg','Ell_1','Ell_inf','Nlp'}
        % reg
        settings_homotopy = settings.homotopy_solver_settings;
        if strcmp(settings.lpec_solver,'Reg')
            % if the lpec solver is scholtes reg, then the direct homotopy is meant
            settings_homotopy.homotopy_parameter_steering = 'Direct';
        elseif strcmp(settings.lpec_solver,'Ell_1')
            settings_homotopy.homotopy_parameter_steering = 'Ell_1';
        elseif strcmp(settings.lpec_solver,'Nlp')
            settings_homotopy.homotopy_parameter_steering = 'None';
        else
            settings_homotopy.homotopy_parameter_steering = 'Ell_inf';
        end
        try
            [result_homotopy,stats] = mpec_homotopy_solver(lpec_homotopy,solver_initalization,settings_homotopy);
        catch
            % keyboard;
            stats.success = false;
            stats.cpu_time_total = 0;
            stats.return_status = 'regularization based lpec solver faild - probbaly infeasible';
        end
        if stats.success
            results.d_lpec = result_homotopy.x;
            results.y_lpec = result_homotopy.x1 >= result_homotopy.x2;
            % results.y_lpec = results.d_lpec(lpec.dims.ind_x1)>=results.d_lpec(lpec.dims.ind_x2);\
            % results.y_lpec = lpec.x_lin(lpec.dims.ind_x1)+results.d_lpec(lpec.dims.ind_x1)>=lpec.x_lin(lpec.dims.ind_x2)+results.d_lpec(lpec.dims.ind_x2);
            results.f_opt = result_homotopy.f;
            stats.lpec_solution_exists = true;
            stats.optimal_solution_found = true;
            
        else
            results.d_lpec = result_homotopy.x;
            results.y_lpec = results.d_lpec(lpec.dims.ind_x1)>=results.d_lpec(lpec.dims.ind_x2);
            results.f_opt = result_homotopy.f;

            stats.lpec_solution_exists = false;
        end
        stats.cpu_time = stats.cpu_time_total;
        stats.solver_message = stats.return_status;
        % stats.nodecount = 0;
    case {"Projected_Gradient"}
        d_lpec = zeros(lpec.dims.n_primal,1);
        d0 = zeros(lpec.dims.n_primal-2*lpec.dims.n_comp,1);
        d1 = zeros(lpec.dims.n_comp,1);
        d2 = d1;
        t_start_pg = tic;
        try
            % solution for ind_x0
            d0(lpec.f(lpec.dims.ind_x0)<0) = min(lpec.ub(lpec.f(lpec.dims.ind_x0)<0)-x_lin_0(lpec.f(lpec.dims.ind_x0)<0),lpec.rho_TR);
            d0(lpec.f(lpec.dims.ind_x0)>0) = max(lpec.lb(lpec.f(lpec.dims.ind_x0)>0)-x_lin_0(lpec.f(lpec.dims.ind_x0)>0),-lpec.rho_TR);
            d0(abs(lpec.f(lpec.dims.ind_x0))<1e-12) = 0;
            for ii = 1:lpec.dims.n_auxiliary
                f_ii = [lpec.f(lpec.dims.ind_x1(ii)); lpec.f(lpec.dims.ind_x2(ii))]; % objectiv of i-th 2d lpec;

                if lpec.rho_TR >= x_lin_1(ii) && x_lin_1(ii) >= 0 && x_lin_2(ii) == 0
                    d_candidates = [0, lpec.rho_TR, -x_lin_1(ii), -x_lin_1(ii);...
                        0, 0, 0, lpec.rho_TR];
                    % case A
                elseif x_lin_1(ii) == 0 && lpec.rho_TR >= x_lin_2(ii) && x_lin_2(ii) >= 0
                    % case B
                    d_candidates = [0, lpec.rho_TR, 0, 0;...
                        0, -x_lin_2(ii), -x_lin_2(ii), lpec.rho_TR];
                elseif x_lin_1(ii)>lpec.rho_TR && x_lin_2(ii) == 0
                    % case C
                    d_candidates = [0, -lpec.rho_TR, lpec.rho_TR;...
                        0, 0,0];
                else
                    % case D
                    d_candidates = [0, 0, 0;...
                        0, -lpec.rho_TR,lpec.rho_TR];
                end
                [val_min,ind_min] = min(d_candidates'*f_ii);
                if val_min == 0
                    d1(ii)  = 0;
                    d2(ii)  = 0;
                else
                    d1(ii) = d_candidates(1,ind_min);
                    d2(ii) = d_candidates(2,ind_min);
                end
            end
            % solution for ind_x1 and ind_x2
            % d_lpec  = [d0;d1;d2];
            d_lpec(lpec.dims.ind_x0) = d0;
            d_lpec(lpec.dims.ind_x1) = d1;
            d_lpec(lpec.dims.ind_x2) = d2;
            f_opt = lpec.f'*d_lpec;
            % y_lpec = d_lpec(lpec.dims.ind_x1)>=d_lpec(lpec.dims.ind_x2);
            % y_lpec = d_lpec(lpec.dims.ind_x1)>=d_lpec(lpec.dims.ind_x2);
            y_lpec = lpec.x_lin(lpec.dims.ind_x1)+d_lpec(lpec.dims.ind_x1)>lpec.x_lin(lpec.dims.ind_x2)+d_lpec(lpec.dims.ind_x2);
            lpec_solution_exists  = true;
            solver_message = 'OPTIMAL';
            stats.optimal_solution_found = true;
        catch
            f_opt = nan;
            y_lpec = nan;
            lpec_solution_exists  = false;
            solver_message = 'INFEASIBLE';
        end
        cpu_total = toc(t_start_pg); % todo add also projection step time;
        results.d_lpec = d_lpec;
        results.y_lpec = y_lpec;
        results.f_opt = f_opt;
        stats.lpec_solution_exists = lpec_solution_exists;
        stats.nodecount = 1;
        stats.cpu_time = cpu_total;
        stats.solver_message = solver_message;
    case "LCLP"
        t_start_lclp = tic;
        [x_opt, stats] = box_lpcc(lb, ub, lpec.f, L, R, lpec.A_eq, lpec.b_eq, lpec.A_ineq, lpec.b_ineq);
        cpu_total = toc(t_start_lclp);
        results.y_lpec = lpec.x_lin(lpec.dims.ind_x1)+x_opt(lpec.dims.ind_x1)>lpec.x_lin(lpec.dims.ind_x2)+x_opt(lpec.dims.ind_x2);
        results.f_opt = stats.f_opt;
        results.d_lpec = x_opt;
        stats.cpu_time = cpu_total;
        stats.lpec_solution_exists = stats.success;
        stats.nodecount = 1;
        if ~stats.success
            stats.solver_message = "FAILED";
        else
            stats.solver_message = "OPTIMAL";
        end
end
end
