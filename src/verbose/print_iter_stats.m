function [] = print_iter_stats(n_major,n_minor,f,h_std,h_comp,solver_type,solver_iters,solver_status,rho_TR,d_norm,cpu_time,accept_trail_step)

max_char = 15;
if length(solver_status) >= max_char
    solver_status  = solver_status(1:max_char);
else
    solver_status  = [repmat(' ',1,1+floor((max_char-length(solver_status))/2)),solver_status,repmat(' ',1, ceil((max_char-length(solver_status))/2)+3)];
end

if ~ischar(n_major)
    n_major = num2str(n_major);
end
% 
try
    if ~ischar(accept_trail_step)
        accept_trail_step = num2str(accept_trail_step);
    end
catch
end
fprintf('|%-7s|%-7d|%-10s|%-10.2e|%-10.2e|%-10.2e|%-10d|%-20s|%-10.2e|%-10.2e|%-5.2e|%-7s|\n', ...
         n_major,n_minor,solver_type,f,h_std,h_comp,solver_iters,solver_status,rho_TR,d_norm,cpu_time,accept_trail_step);
end