function active_set_guess = find_active_sets_piece_nlp(x_trail,nabla_f,y_lpec,dims,settings,tol)
active_set_guess = struct;
% Note that this active set selection assume that the complementarities are feasible.

% get the basic active sets as in the piece nlps
% tol = tol*10;
% addaptive tolerance;
% tol = 1.1*max(min([abs(x_trail(dims.ind_x1)),abs(x_trail(dims.ind_x2))],[],2))+tol;
% active_set_guess.I_00 = 0;
% active_set_guess.I_0_plus = 0;
% active_set_guess.I_plus_0 = 0;
% ii_active_set = 1;
% ii_max = 10;
% tol = settings.tol_active;
if dims.n_comp > 0
    % while sum(active_set_guess.I_00+active_set_guess.I_0_plus+active_set_guess.I_plus_0) < dims.n_comp && ii_active_set <= ii_max
        active_set_guess.ind_x1_active = abs(x_trail(dims.ind_x1)) <= tol; % A1(x) in Kirches2022
        active_set_guess.ind_x2_active = abs(x_trail(dims.ind_x2)) <= tol; % A2(x) in Kirches2022
        active_set_guess.ind_biactive = active_set_guess.ind_x1_active & active_set_guess.ind_x2_active; % D(x)  = A1(X) \cap A2(x) in Kirches2022
        active_set_guess.ind_x1_strictly_active = x_trail(dims.ind_x1) > tol & x_trail(dims.ind_x2) <= tol; % A1+(x) in Kirches2022
        active_set_guess.ind_x2_strictly_active = x_trail(dims.ind_x2) > tol & x_trail(dims.ind_x1) <= tol;  % A2+(x) in Kirches2022
        % diff notation
        active_set_guess.I_plus_0 = active_set_guess.ind_x1_strictly_active;
        active_set_guess.I_0_plus = active_set_guess.ind_x2_strictly_active;
        active_set_guess.I_00 = active_set_guess.ind_biactive;
    %     ii_active_set  = ii_active_set + 1;
    %     tol = 5*tol;
    %     % if ii_active_set  > 2
    %     %     keyboard
    %     % end
    % end
end

switch settings.piece_nlp_strategy
    case 'TNLP'
        % take the base as it is from the base;
        if sum(active_set_guess.I_00+active_set_guess.I_0_plus+active_set_guess.I_plus_0) < dims.n_comp
            if settings.debug_mode_on
                % keyboard;
            end
            active_set_guess = addaptive_active_set_identification(x_trail,dims,tol);
            % warning('Active set not extract from lpec solution - fallback: solving a BNLP.')
            % active_set_guess.I_plus_0 = y_lpec==1;
            % active_set_guess.I_0_plus = y_lpec==0;
            % active_set_guess.I_00 = [];
        end
    case 'BNLP_integer'
        active_set_guess.I_plus_0 = y_lpec==1;
        active_set_guess.I_0_plus = y_lpec==0;
        active_set_guess.I_00 = [];
    case "Gradient1"
        % split I_00
        active_set_guess.I_0_plus(active_set_guess.I_00) = nabla_f(dims.ind_x1(active_set_guess.I_00))>=nabla_f(dims.ind_x2(active_set_guess.I_00));
        active_set_guess.I_plus_0(active_set_guess.I_00) = nabla_f(dims.ind_x1(active_set_guess.I_00))<nabla_f(dims.ind_x2(active_set_guess.I_00));
    case "Gradient2"
        % split I_00
        active_set_guess.I_0_plus(active_set_guess.I_00) = nabla_f(dims.ind_x1(active_set_guess.I_00))>nabla_f(dims.ind_x2(active_set_guess.I_00));
        active_set_guess.I_plus_0(active_set_guess.I_00) = nabla_f(dims.ind_x1(active_set_guess.I_00))<=nabla_f(dims.ind_x2(active_set_guess.I_00));

end

end