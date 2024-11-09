function [] = print_iter_header_objective()
line_str = repmat('-',1,124);
line_str = ['\n' line_str];
fprintf(line_str)
fprintf('\n|%-7s|%-7s|%-9s|%-9s|%-10s|%-15s|%-15s|%-15s|%-10s|%-10s|%-10s|%-10s|%-10s\n',...
           'M.iter', 'm.iter', 'f(x^k)', 'f(x^nlp)', 'Delta f', '||nabla_f||' , 'nabla_f''d_nlp', 'nabla_f''d_lpec', '||d_nlp||', '||d_lpec||', 'rho_TR', 'NLP resolve');
fprintf(line_str)
end

