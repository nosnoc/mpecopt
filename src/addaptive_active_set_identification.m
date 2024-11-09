function active_set_guess = addaptive_active_set_identification(x_trail,dims,tol)
active_set_guess = struct;
% tol = tol*10;
ii_active_set  = 1;
max_iter = 10;
active_set_guess.I_plus_0 = 0;
active_set_guess.I_0_plus = 0;
active_set_guess.I_00 = 0;

if dims.n_comp > 0
    while sum(active_set_guess.I_00+active_set_guess.I_0_plus+active_set_guess.I_plus_0) < dims.n_comp && ii_active_set <= max_iter
        active_set_guess.ind_x1_active = abs(x_trail(dims.ind_x1)) <= tol; % A1(x) in Kirches2022
        active_set_guess.ind_x2_active = abs(x_trail(dims.ind_x2)) <= tol; % A2(x) in Kirches2022
        active_set_guess.ind_biactive = active_set_guess.ind_x1_active & active_set_guess.ind_x2_active; % D(x)  = A1(X) \cap A2(x) in Kirches2022
        active_set_guess.ind_x1_strictly_active = x_trail(dims.ind_x1) > tol & x_trail(dims.ind_x2) <= tol; % A1+(x) in Kirches2022
        active_set_guess.ind_x2_strictly_active = x_trail(dims.ind_x2) > tol & x_trail(dims.ind_x1) <= tol;  % A2+(x) in Kirches2022
        % diff notation
        active_set_guess.I_plus_0 = active_set_guess.ind_x1_strictly_active;
        active_set_guess.I_0_plus = active_set_guess.ind_x2_strictly_active;
        active_set_guess.I_00 = active_set_guess.ind_biactive;
        ii_active_set  = ii_active_set + 1;
        tol = 5*tol;
    end
end
active_set_guess.tol = tol;

end