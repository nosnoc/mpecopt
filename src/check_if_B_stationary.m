function [success,stopping_criterion_fullfiled,solver_message] = check_if_B_stationary(d_lpec,f_lin_opt,nabla_f,h_total,success,stopping_criterion_fullfiled,solver_message,settings)

% check if feasible and if d = 0 is a valid solution
% if objective zero (either if cost gradient zero, or solution leads to it) = then set step to zero => B stationarity
if (h_total < settings.tol) && ((abs(f_lin_opt) <= settings.tol_B_stationarity || norm(nabla_f) < settings.tol_B_stationarity))     
                % if the current point is feasible, and the objective is zero, then d = 0 is also a solution of the lpec (occurs if a solution is not on the verties of the lp)
                d_lpec = d_lpec*0; 
end
% check if satisfies 
   if norm(d_lpec) <= settings.tol_B_stationarity
                stopping_criterion_fullfiled = true;     % B-stationary point found, optimal solution found!
                solver_message = 'B-stationary point found sucessfully.';
                success = true; 
   else
               % status quol
    end
end