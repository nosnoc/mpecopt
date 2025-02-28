classdef TestPhaseTwoOpts < matlab.unittest.TestCase
    properties (TestParameter)
        % Only test the 
        piece_nlp_strategy = {'TNLP', 'BNLP_integer', 'BNLP_Gradient1', 'BNLP_Gradient2'};
    end

    methods (Test, ParameterCombination = 'exhaustive')
        function test_initialization_strategy(tc,piece_nlp_strategy)
            import casadi.*
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tc.applyFixture( ...
                SuppressedWarningsFixture("optim:intlinprog:IgnoreX0"))

            x1 = SX.sym('x1');
            x2 = SX.sym('x2');
            x3 = SX.sym('x3');
            x4 = SX.sym('x4');
            x5 = SX.sym('x5');
            x6 = SX.sym('x6');
            x7 = SX.sym('x7');
            x8 = SX.sym('x8');

            w = [x1;x2;x3;x4;x5;x6;x7;x8];

            p = SX.sym('p');
            f = (x1-5)^2+(2*x2+1)^2;
            g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
                3*x1-x2-3-x6;...
                -x1+0.5*x2+4-x7;...
                -x1-x2+7-x8];
            G = [x6;x7;x8];
            H = [x3;x4;x5];
            x0 = [0;0;1;0;0;0;1;1];
            lbw = zeros(8,1);
            ubw = inf*ones(8,1);

            lbg = zeros(4,1);
            ubg = zeros(4,1);

            mpec = struct('x', w, 'f', f, 'g', g,'p',p,'G',G ,'H',H);
            solver_initalization = struct('x0', x0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg, 'ubg',ubg,'p0',1, 'y0', [0;1;1]);

            opts = mpecopt.Options();
            opts.piece_nlp_strategy = piece_nlp_strategy;
            %opts.verbose_solver = false;
            %opts.verbose_summary = false;

            solver = mpecopt.Solver(mpec, opts);
            [sol,stats] = solver.solve(solver_initalization);
        end
    end
end
