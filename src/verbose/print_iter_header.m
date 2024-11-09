function [] = print_iter_header()

line_str = repmat('-',1,132);
fprintf(['\n' line_str])
fprintf('\n|%-7s|%-7s|%-10s|%-10s|%-10s|%-10s|%-10s|%-20s|%-10s|%-10s|%-7s|%-7s|\n', ...
         'M.iter', 'm.iter', 'Subprob', 'objective', 'h_std', 'h_comp', 'Iters', 'Status','rho TR/tau','||d||/D_f','Time (s)','Step Acc.');
fprintf([line_str '\n'])

end

