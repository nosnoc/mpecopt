function [] = print_iter_summary_minlp(f,h_std,h_comp,solver_message,multiplier_based_stationarity,b_stationarity,n_biactive,f_lpec,rho_TR)

print_iter_line()
fprintf('\n');
fprintf('MPEC MINLP reformulation optimizer:\t %s\n',solver_message);
fprintf('Objective....................:\t %2.6e\n',f);
fprintf('Std. constraint violation....:\t %2.6e\n',h_std);
fprintf('Complementarity residual.....:\t %2.6e\n',h_comp);
% fprintf('Solver message...............:\t %s\n',solver_message)
fprintf('Mult. based stationarity.....:\t %s\n',multiplier_based_stationarity);
fprintf('B stationarity...............:\t %s\n',mat2str(b_stationarity));
fprintf('Biactive constraints.........:\t %d\n',n_biactive);
fprintf('nabla_f(x)^T d...............:\t %d\n',f_lpec);
fprintf('Final rho_TR.................:\t %d\n',rho_TR);
print_iter_line()
fprintf('\n');
end