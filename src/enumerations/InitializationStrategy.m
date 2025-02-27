classdef InitializationStrategy
    enumeration
        RelaxAndProject  % solve a scholtes relaxed problem, then an LPEC to get an active set, solve BNLP after and take this as x0;
        TakeInitialGuessDirectly % use the provided initial guess and just start with phase II (make sense if x0 is a feasible point of some BNLP)
        TakeInitialGuessActiveSet % take the BNLP generated from the active set at x0;
        TakeProvidedActiveSet % If a binary y encoding is I_1 and I_2 is avilable, then take this BNLP, otherwise do same as TakeInitialGuessActiveSet;
        RandomActiveSet % just take the user provided x0;
        AllBiactive % set all complementarities x_1 = 0, x_2 = 0 and solve this bnlp;
        FeasibilityEll1General % only ell1 penalization of  g_lb <= g(x) <= g_ub; G(x) = x1 and H(x) = x2; Note that equality constraint can be relaxed with both one and two slacks.
        FeasibilityEllInfGeneral % only ell_infity penalization of  g_lb <= g(x) <= g_ub G(x) = and x1 and H(x) = x2;
        % FeasibilityEllInfAll, % only ell1 penalization of  g_lb <= g(x) <= g_ub, G(x) = x1 and H(x) = x2;  and   x0_lb <= x0 <= x0_ub
        % FeasibilityEll1All, % only ell1 penalization of  g_lb <= g(x) <= g_ub, G(x) = x1 and H(x) = x2;  and   x0_lb <= x0 <= x0_ub
    end
    
end

