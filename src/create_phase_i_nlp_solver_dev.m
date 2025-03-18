function   [solver_relaxed,solver_initialization_relaxed] = create_phase_i_nlp_solver_dev(mpec_casadi,solver_initialization,settings,dims)

% Remark: This could also be done directly with the scholtes solver, however I do it explicitly to avoid the preprocess overhead.
import casadi.*
% symbolics
sigma_relaxed = SX.sym('sigma_relaxed',1); 
x_relaxed = mpec_casadi.x; 
x1 = mpec_casadi.x1;
x2 = mpec_casadi.x2;
f_relaxed = mpec_casadi.f;
g = mpec_casadi.g;
x_k_relaxed = solver_initialization.x0;
% Parameters
p_relaxed = [mpec_casadi.p;sigma_relaxed];  
p0_relaxed = [solver_initialization.p0; settings.relax_and_project_sigma0];

% relaxed complementarity constraints;
g_comp_relaxed = []; 
lbg_comp_relaxed = [];  
ubg_comp_relaxed = [];

% bounds
lbx_relaxed = solver_initialization.lbx; 
ubx_relaxed = solver_initialization.ubx;
lbg_relaxed = solver_initialization.lbg; 
ubg_relaxed = solver_initialization.ubg;

switch settings.relax_and_project_homotopy_parameter_steering
    case "Direct"
        if settings.relax_and_project_comps_aggregated
            g_comp_relaxed = x1'*x2-sigma_relaxed*dims.n_comp;
            lbg_comp_relaxed = -inf;
            ubg_comp_relaxed = 0;
        else
            g_comp_relaxed = x1.*x2-sigma_relaxed;
            lbg_comp_relaxed = -inf*ones(dims.n_comp,1);
            ubg_comp_relaxed = 0*ones(dims.n_comp,1);
        end
    case "Ell_1"
        f_relaxed = f_relaxed+(x1'*x2)*(sigma_relaxed)^(-1);
    case "Ell_inf"
        % x_k_relaxed = project_to_bounds(x_k_relaxed ,lbx,ubx,dims);
        s_eleastic = SX.sym('s_eleastic',1);
        f_relaxed = f_relaxed+(s_eleastic)*(sigma_relaxed)^(-1);
        x_relaxed = [x_relaxed;s_eleastic];
        lbx_relaxed = [lbx_relaxed;0];
        ubx_relaxed = [ubx_relaxed;max(10,settings.relax_and_project_sigma0*10)];
        x_k_relaxed = [x_k_relaxed;settings.relax_and_project_sigma0];
        if settings.relax_and_project_comps_aggregated
            g_comp_relaxed = x1'*x2-s_eleastic*dims.n_comp;
            lbg_comp_relaxed = -inf;
            ubg_comp_relaxed = 0;
        else
            g_comp_relaxed = x1.*x2-s_eleastic;
            lbg_comp_relaxed = -inf*ones(dims.n_comp,1);
            ubg_comp_relaxed = 0*ones(dims.n_comp,1);
        end
end

% ind_comp = length(g)+1;
g_relaxed = [g;g_comp_relaxed]; 
lbg_relaxed = [lbg_relaxed;lbg_comp_relaxed]; 
ubg_relaxed = [ubg_relaxed;ubg_comp_relaxed];
% ind_comp = ind_comp:1:length(g_relaxed);
% Output solver initalization
solver_initialization_relaxed.x0 = x_k_relaxed;
solver_initialization_relaxed.lbx = lbx_relaxed;
solver_initialization_relaxed.ubx = ubx_relaxed;
solver_initialization_relaxed.lbg = lbg_relaxed;
solver_initialization_relaxed.ubg = ubg_relaxed;
solver_initialization_relaxed.p0 = p0_relaxed;

%% --------- create solver for Phase I -----------------------------------
nlp_relaxed = struct('x', x_relaxed,'f', f_relaxed,'g', g_relaxed,'p',p_relaxed);
solver_relaxed = nlpsol('solver_relaxed', 'ipopt', nlp_relaxed, settings.settings_casadi_nlp);
end