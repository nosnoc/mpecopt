function [] = print_iter_stats_objective(k,l_k,f_k,f_k_trail,nabla_f,inner_nonlin,inner_lin,d_nlp_k,d_lpec_k,rho_TR_k_l,resolve_nlp)
    fprintf('\n|%-7d|%-7d|%-9.2e|%-9.2e|%-10.2e|%-15.2e|%-15.2e|%-15.2e|%-10.2e|%-10.2e|%-10.2e|%-10d',...
                 k,  l_k,f_k,f_k_trail,f_k-f_k_trail,norm(nabla_f),inner_nonlin,inner_lin,norm(d_nlp_k),norm(d_lpec_k),rho_TR_k_l,resolve_nlp);

% inner_nonlin = nabla_f'*d_nlp_k;
% inner_lin = nabla_f'*d_lpec_k;
end