function [] = print_iter_summary(f,h_std,h_comp,solver_message,multiplier_based_stationarity,b_stationariry,n_biactive,f_lpec,rho_TR)

% print_iter_line()
line_str = repmat('-',1,130);
line_str = ['\n' line_str '\n' ];
fprintf(line_str)
fprintf('MPECopt:\t %s\n',solver_message);
fprintf('Objective....................:\t %2.6e\n',f);
fprintf('Std. constraint violation....:\t %2.6e\n',h_std);
fprintf('Complementarity residual.....:\t %2.6e\n',h_comp);
% fprintf('Solver message...............:\t %s\n',solver_message)
fprintf('Mult. based stationarity.....:\t %s\n',multiplier_based_stationarity);
fprintf('B stationarity...............:\t %s\n',mat2str(b_stationariry));
fprintf('Biactive constraints.........:\t %d\n',n_biactive);
fprintf('nabla_f(x)^T d...............:\t %d\n',f_lpec);
fprintf('Final rho_TR.................:\t %d\n',rho_TR);
print_iter_line()
fprintf('\n');
end