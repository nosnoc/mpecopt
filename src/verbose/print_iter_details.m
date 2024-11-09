function [] = print_iter_details(n_nlp,n_lpec,rho_TR_iter,create_solver_time,total_cpu_time,presolve_cpu_time,main_loop_cpu_time)

if isempty(rho_TR_iter)
    rho_TR_iter = nan;
end
line_str = repmat('-',1,130);
line_str = ['\n' line_str '\n' ];
fprintf(line_str)
% TODO separate in preproces and in mailoop data;
fprintf('NLPs solved..................:\t %d\n',n_nlp);
fprintf('LPECs solved.................:\t %d\n',n_lpec);
fprintf('Max. TR Radius...............:\t %2.6e\n',max(rho_TR_iter));
fprintf('Min. TR Radius...............:\t %2.6e\n',min(rho_TR_iter));
% todo create time conversion
fprintf('Create solver time (s).......:\t %2.2f\n',create_solver_time);
fprintf('Presolve time (s)............:\t %2.2f\n',presolve_cpu_time);
fprintf('Main iterations time (s).....:\t %2.2f\n',main_loop_cpu_time);
fprintf('Total time (s)................:\t %2.2f\n',total_cpu_time);

% print_iter_line()

end