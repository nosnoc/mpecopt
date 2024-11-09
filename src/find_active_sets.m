function active_set_guess = find_active_sets(x_trail,dims,tol)
active_set_guess = struct;
% Note that this active set selection assume that the complementarities are feasible.
if dims.n_comp > 0

    active_set_guess.ind_x1_active = abs(x_trail(dims.ind_x1)) <= tol; % A1(x) in Kirches2022
    active_set_guess.ind_x2_active = abs(x_trail(dims.ind_x2)) <= tol; % A2(x) in Kirches2022
    active_set_guess.ind_biactive = active_set_guess.ind_x1_active & active_set_guess.ind_x2_active; % D(x)  = A1(X) \cap A2(x) in Kirches2022
    active_set_guess.ind_x1_strictly_active = x_trail(dims.ind_x1) > tol & x_trail(dims.ind_x2) <= tol; % A1+(x) in Kirches2022
    active_set_guess.ind_x2_strictly_active = x_trail(dims.ind_x2) > tol & x_trail(dims.ind_x1) <= tol;  % A2+(x) in Kirches2022

    active_set_guess.I_plus_0 = active_set_guess.ind_x1_strictly_active;
    active_set_guess.I_0_plus = active_set_guess.ind_x2_strictly_active;
    active_set_guess.I_00 = active_set_guess.ind_biactive;
end

end