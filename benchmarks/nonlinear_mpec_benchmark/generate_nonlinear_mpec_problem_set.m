function mpecs = generate_nonlinear_mpec_problem_set(problem_set_name,settings,dimensions)
% This function generats a set of generic nonlinear MPECs, where are problem variables and function expressions are
% in the CasADi symbolic framework.
% min cx+dy
% s.t.  Ax+By >= f
%       0<= y _|_ q+Nx+My >=0
%       x>=0
%       -Delta <=(x,y)<= Delta
% Reformulate into vertical form
% w = (x,y,z)
% min f(x,y)+cx+dy
% s.t.  Ax+By >= f          [A,B,0]w-f >=0   % (some of these constraint can be copied, to violate LICQ, but keep MFCQ satisfied)
%       q+Nx+My-Iz = 0;     [N,M,-I]w+q = 0 % (some of the constraint can be squared, or to the power of four to make the problem nonlinear)
%       0<= y _|_ z >=0
%       x>=0
%       0 <=w<= w_ub

% n_x number of noncomplementarity variablex (x0)
% n_y number of complementarity variablex (x1,x3)
% n = n_x+2n_y total number of variables

% LPEC generator is inspired by:
%
% [1] Hu, Jing, John E. Mitchell, Jong-Shi Pang, Kristin P. Bennett, and Gautam Kunapuli.
% "On the global solution of linear programs with linear complementarity constraints."
% SIAM Journal on Optimization 19, no. 1 (2008): 445-471.
%
% The problem generation is described nicely on page 22 of :
% [2] Jara-Moroni, Francisco, John E. Mitchell, Jong-Shi Pang, and Andreas Wächter.
% "An enhanced logical benders approach for linear programs with complementarity constraints."
% Journal of Global Optimization 77 (2020): 687-714.

% Nonlinear Objective functions ispired by:
%
% [3] N. Andrei, An unconstrained optimization test functions collection, Advanced Modeling and Optimization,
% 10 (2008), pp. 147–161, https://camo.ici.ro/journal/vol10/v10a10.pdf.

% by Armin Nurkanovic, University of Freiburg, Mai, 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import casadi.*
rng(1,"twister"); % to have reproducable problem sets
unfold_struct(settings,'caller')
unfold_struct(dimensions,'caller')


if ~isfield(settings,'nnz_bounded_by_dim')
    nnz_bounded_by_dim = false;
end

if ~isfield(settings,'n_comp_min')
    n_comp_min = 500; % after this trashold addaptiv bound used
end


if ~isfield(settings,'inv_cond_num')
    % inv_cond_num = 1e-1;
    inv_cond_num = [];
end

if ~isfield(settings,'nnz_factor')
    nnz_factor = 1;
end


if ~isfield(settings,'casadi_variable_type')
    casadi_variable_type = 'SX'; % default;
else
    if ~ismember(settings.casadi_variable_type,{'SX','MX'})
        error('casadi_variable_type must be = ''SX'' or ''MX'' ') 
    end
end

eps_prec = 10^(-settings.n_digits);

if  settings.round_all_data
    n_digits_data = settings.n_digits;
else
    n_digits_data = 16;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_objectives = length(objective_functions);
if settings.random_problem_sizes
    % problems size across all
    latexify_plot();
    N_problems = N_objectives*N_rand_prob;
    n_x_max = max(N_problems+n_x_min,n_x_max);
    % n_x_vec = n_x_min+randperm(n_x_max-n_x_min,N_problems); % num of non comp vars
    % Create a weighted distribution favoring smaller values
    weights = linspace(5, 1, n_x_max-n_x_min+1);  % Higher weights for smaller indices
    weights = weights / sum(weights);  % Normalize
    % Generate biased random indices
    biased_indices = randsample(1:(n_x_max-n_x_min+1), N_problems, true, weights);
    n_x_vec = n_x_min + biased_indices - 1;
    n_y_vec = round(n_x_vec/n_fraction_of_x); % num of comp vars
    rand_ineq_scale = settings.n_ineq_lb+(settings.n_ineq_ub-settings.n_ineq_lb)*(rand(1,N_problems));
    n_ineq_vec = round(rand_ineq_scale.*n_x_vec);
    % n_ineq_vec = round((settings.n_ineq_lb+(settings.n_ineq_ub-settings.n_ineq_lb)*(rand(1,N_problems)).*n_x_vec))*3;

    N_per_problem = N_rand_prob;
    if 1
    figure
    n_var = n_x_vec+2*n_y_vec; % total num of vars, non comp, comp1 + comp2 + lift
    scatter(n_var,n_ineq_vec*(1+settings.s_ineq_copy));
    hold on
    scatter(n_var,n_x_vec);
    hold on
    scatter(n_var,n_y_vec);
    xlabel('$n+2m$ - number of variables')
    ylabel('Number of constraints')
    legend({'Inequality constraints','Equality constraints', 'Comp. constraints'},'Location','northwest')
    axis equal
    end
    % subplot(122)
    % scatter(n_var,n_y_vec);
    % hold on
    % scatter(n_var,n_ineq_vec*(1+settings.s_ineq_copy)+n_x_vec);
    % scatter(n_var,n_x_vec);
    % xlabel('$n+2m$ - number of variables')
    % legend({'$m$ number of complementairt variables','Total num. of constraints', 'Non. comp. vars.'},'Location','northwest')
    % axis equal
    % elseif settings.random_problem_sizes_individual
    %     n_x_vec = round(n_x_min+n_x_max*rand(N_rand_prob,1));
    %     n_x_vec = sort(n_x_vec);
    %     n_y_vec = round(n_x_vec/n_fraction_of_x);
    %     n_ineq_vec = round((settings.n_ineq_lb+(settings.n_ineq_ub-settings.n_ineq_lb)*(rand(1,N_problems)).*n_x_vec));
    %
    %     n_x_vec = repmat(n_x_vec',1,N_objectives);
    %     n_y_vec = repmat(n_y_vec',1,N_objectives);
    %     n_ineq_vec = repmat(n_ineq_vec',1,N_objectives);
    %
    %     % Vectors stored in dimensions, example for the format;
    %     % n_x_vec = [20,50,100];
    %     % n_y_vec = [20,50,100];
    %     % k_vec = [10,35,90];
    %     N_per_problem = N_rand_prob;
else
    N_per_problem = length(n_x_vec);
    n_x_vec = repmat(n_x_vec',1,N_objectives);
    n_y_vec = repmat(n_y_vec',1,N_objectives);
    n_ineq_vec = repmat(n_ineq_vec',1,N_objectives);
end

fprintf('Largerst problem: n_var = %d, n_comp = %d \n',max(n_x_vec+2*n_y_vec), max(n_y_vec))
n_comp_max = max(n_y_vec);

% range for feasible points
range_x = [0,1];
range_y = [0,1];

% range for objective gradient
range_c = [0, 1]./rescale_factor;
range_d = [1, 3]./rescale_factor;

% range for problem matrices
if ~isfield(settings,'range_s_density')
    range_s_density = [0.01 0.1];
else
    range_s_density = settings.range_s_density;
end


if ~isfield(settings,'adaptive_density_bounds')
    adaptive_density_bounds = 0;
else
    adaptive_density_bounds = settings.adaptive_density_bounds;
end

range_A = [0,2]./rescale_factor; % uniform
range_B = [0,2]./rescale_factor; % uniform
% range_A = [-1,1]./rescale_factor; % uniform
% range_B = [-1,1]./rescale_factor; % uniform
range_d1 = [0,2]./rescale_factor; % uniform
range_d2 = [0,2]./rescale_factor; % uniform
range_N = [-1,1]./rescale_factor; % uniform
range_E = [-1,1]./rescale_factor; % uniform
range_E = [3,4]./rescale_factor; % uniform
range_q = [-20,10]./rescale_factor; % uniform
range_grad = [-10,10];
range_ubw = [1e1,0.5e3];

mpecs = [];
n_eq_vec = [];
% n_ineq_vec = [];
n_var_vec = [];
n_comp_vec = [];


% Generate values from the uniform distribution on the interval (a, b). r = a + (b-a).*rand(100,1);
for kk = 1:length(objective_functions)
    objective_type = objective_functions{kk};
    for jj = 1:N_per_problem
        % dimensions
        % ii_jj = N_per_problem*(jj-1)+ii;
        ii_kk = N_per_problem*(kk-1)+jj;
        n_x = n_x_vec(ii_kk);
        n_y = n_y_vec(ii_kk);
        n_ineq = n_ineq_vec(ii_kk);
        % Ax+By >= f;  x is of n_x, y is of n_y, f is of n_ineq;
        n = n_x+n_y+n_y; % total number of vairables  [x,y,z] ; z = q+Nx+My;
        n_comp = n_y;
        n_non_lifted = n_x+n_y;
        % store dimensions:
        n_var_vec = [n_var_vec, n];
        n_comp_vec = [n_comp_vec, n_comp];
        % Define symbolic variables:
        if isequal(casadi_variable_type,'SX')
            x = SX.sym('x', n_x);
            y = SX.sym('y', n_y);
            z = SX.sym('z', n_y);
        else
            x = MX.sym('x', n_x);
            y = MX.sym('y', n_y);
            z = MX.sym('z', n_y);
        end
        w = [x;y;z];
        v = [x;y];

        if settings.include_lifted_variables_in_objective
            n_obj = n_x+n_y+n_y;
        else
            n_obj = n_non_lifted;
        end

        c = range_c(1)+(range_c(2)-range_c(1)).*rand(n_x,1);
        d = range_d(1)+(range_d(2)-range_d(1)).*rand(n_y,1);


        f_lin = [c;d]'*v;
        f = alpha_lin*f_lin;

        if settings.variable_density
            % Addaptive bounds on the densitiy, as for very larg problems less sparse problems can take hours to create an expression.
            % graph
            if adaptive_density_bounds
                % for < 1000 keep fixed, for larger shrink linearly lower and upper bounds
                % range_s_density_lb = max(0.0025, min(range_s_density(1), range_s_density(1) - (n_comp-n_comp_min)*(0.0075/(n_comp_max-n_comp_min))));
                % range_s_density_ub = max(0.005,  min(range_s_density(2), range_s_density(2) - (n_comp-n_comp_min)*(0.045/(n_comp_max-n_comp_min))));
                if n_comp <= n_comp_min
                    % Below minimum: use original range_s_density bounds
                    range_s_density_lb = range_s_density(1);
                    range_s_density_ub = range_s_density(2);
                elseif n_comp >= n_comp_max
                    % At or above maximum: use target bounds
                    range_s_density_lb = 0.0025;
                    range_s_density_ub = 0.005;
                else
                    % Between min and max: linear interpolation
                    % Calculate interpolation factor (0 at n_comp_min, 1 at n_comp_max)
                    alpha = (n_comp - n_comp_min) / (n_comp_max - n_comp_min);

                    % Interpolate: shrink from original wide bounds to tighter target bounds
                    range_s_density_lb = range_s_density(1) + alpha * (0.0025 - range_s_density(1));
                    range_s_density_ub = range_s_density(2) + alpha * (0.005 - range_s_density(2));
                end
                s_density_A_B = range_s_density_lb +(range_s_density_ub-range_s_density_lb).*rand(1);
            else
                range_s_density_lb = range_s_density(1);
                range_s_density_ub = range_s_density(2);
                s_density_A_B = range_s_density_lb +(range_s_density_ub-range_s_density_lb).*rand(1);
            end
        end

        switch objective_type
            case 'Linear'
                % Objective
                % c = range_c(1)+(range_c(2)-range_c(1)).*rand(n_x,1);
                % d = range_x(1)+(range_x(2)-range_x(1)).*rand(n_y,1);
                % f = f+[c;d]'*v;

            case 'Quadratic_psd'
                % Quadratic
                % ... generate sparse Hessian matrix H (and make sure it it psd if needed)
                % Quadratic
                if nnz_bounded_by_dim 
                    if adaptive_density_bounds
                        H  = full( sprandsym(n_obj, s_density_A_B) );     % for smaller matrices use provided density bounds, for larger use sparse
                    else
                        if isempty(inv_cond_num)
                            H  = full( sprandsym(n_obj, n_obj/n_obj.^2) );               % ... full matrix for easier printing
                        else
                            H  = full( sprandsym(n_obj, n_obj/n_obj.^2,inv_cond_num) );               % ... full matrix for easier printing
                        end
                    end
                else
                    H  = full( sprandsym(n_obj, s_density_A_B) );               % ... full matrix for easier printing    
                end
                H  = round( H , n_digits_data);                        % ... round matrix to nDigits
                min_eig = max( 0, -min( eig( H ) - eps_prec) );  % ... -epsPrec to ensure H is p.s.d.
                H = H + 2*round( min_eig*eye(n_obj), n_digits);     % ... make trailing matrix p.s.d.
                grad_vec  = round( range_grad(1) + rand(n_obj,1) * (range_grad(2) - range_grad(1)) , n_digits );
                f = 0.5*w(1:n_obj)'*H*w(1:n_obj)+grad_vec'*w(1:n_obj);
                % f = 0.5*v'*H*v+grad_vec'*w(1:n_obj); % +[c;d]'*v;
            case 'Quadratic_ind'
                if settings.nnz_bounded_by_dim
                    if isempty(inv_cond_num)
                        H  = full( sprandsym(n_obj, n_obj/n_obj.^2) );               % ... full matrix for easier printing
                    else
                        H  = full( sprandsym(n_obj, n_obj/n_obj.^2,inv_cond_num) );               % ... full matrix for easier printing
                    end
                else
                    H  = full( sprandsym(n_obj, s_density_A_B) );               % ... full matrix for easier printing    
                end
                H  = round( H , n_digits_data);
                grad_vec  = round( range_grad(1) + rand(n_obj,1) * (range_grad(2) - range_grad(1)) , n_digits );
                f = 0.5*w(1:n_obj)'*H*w(1:n_obj)+grad_vec'*w(1:n_obj);
                % f = 0.5*w(1:n_obj)'*H*w(1:n_obj)+[c;d]'*v;
                % f = 0.5*v'*H*v+grad_vec'*v; % +[c;d]'*v;
            case 'Fletcher'
                for i = 1:n_obj-1
                    f = f+(w(i+1)-w(i)+1-w(i)^2)^2;
                end
            case 'Himmelblau'
                for i = floor(n_obj/2)
                    f = f+(w(2*i-1)-w(2*i)-11)^2+(w(2*i-1)+w(2*i)^2-7)^2;
                end

            case 'McCormick'
                for i = 1:n_obj-1
                    f = f+(-1.5*w(i)+2.5*w(i+1)+1+(w(i)-w(i+1))^2+sin(w(i)+w(i+1)));
                end
            case 'Powell'
                for i = 1:floor(n_obj/4)
                    f = f+ ( (w(4*i-3)+10*w(4*i-2))^2 + 5*(w(4*i-1)-w(4*i))^2 + (w(4*i-2)-2*w(4*i-1))^4 + 10*(w(4*i-3)-w(4*i))^4) ;
                end
            case 'Trigonometric'
                % for i = 1:n_obj-1
                %     f = f+(-1.5*cos(w(i))+2.5*tanh(w(i+1)/0.1)+1+(w(i)-w(i+1))^4+exp(w(i)+w(i+1)));
                % end
                for i = 1:n_obj
                    f = f+ ((n_obj-sum(cos(w))) + i*(1-cos(w(i)))-sin(w(i)))^2;
                end
            case 'Rosenbrock'
                for i = 1:n_obj-1
                    f = f+100*(w(i+1)-w(i)^2)^2+(1-w(i))^2;
                end
            case 'Raydan1'
                for i = 1:n_obj
                    f = f+i/10*(exp(w(i))-w(i));
                end
            case 'Raydan2'
                for i = 1:n_obj
                    f = f+(exp(w(i))-w(i));
                end

            case 'Diagonal3'
                for i = 1:n_obj
                    f = f+(exp(w(i))-i*sin(w(i)))/1e3;
                end
            case 'Diagonal4'
                for i = 1:floor(n_obj/2)
                    f = f+0.5*(w(2*i-1)^2+100*w(2*i)^2);
                end
            case 'Diagonal5'
                for i = 1:n_obj
                    f = f+log(exp(w(i))-exp(-w(i)));
                end
            case 'Extended_Trigiaonal'
                for i = 1:floor(n_obj/2)
                    f = f+(w(2*i-1)+w(2*i)-3)^2+(w(2*i-1)-w(2*i)-3)^4;
                end
            case 'Three_Exponential_Terms'
                for i = 1:floor(n_obj/2)
                    f = f+exp(w(2*i-1)+3*w(2*i)-0.1)+exp(w(2*i-1)-3*w(2*i)-0.1)+exp(-w(2*i-1)-0.1);
                end
            case 'Generalized_PSC1'
                for i = 1:n_obj-1
                    f = f+(w(i)^2+w(i+1)^2+w(i)*w(i+1))^2+sin(w(i))^2+cos(w(i))^2;
                end
            case 'Extended_PSC1'
                for i = 1:floor(n_obj/2)
                    f = f+(w(2*i-1)^2+w(2*i)^2+w(2*i-1)*w(2*i))^2+sin(w(2*i-1))^2+cos(w(2*i))^2;
                end
            case 'Fletcvb3'
                p = 1/1e8;
                c = 1;
                f = f+ 0.5*p*(w(1)^2+w(end)^2);
                for i = 1:n_obj-1
                    f = f+0.5*p*(w(i)-w(i+1))^2;
                end
                for i = 1:n_obj
                    h = 1/(i+1);
                    f = f-(((p*(h^2+2))/h^2)*w(i)+(c*p/h^2)*cos(w(i)));
                end
            case 'Bdqrtic'
                for i = 1:n_obj-4
                    f = f+ (-4*w(i)+3)^2+(w(i)^2+2*w(i+1)^2+3*w(i+2)^2+4*w(i+3)^2+5*w(i+4)^2);
                end
            case 'Tridia'
                f = f+(w(1)-1)^2;
                for i = 2:n_obj
                    f = f+ i*(2*w(i)-w(i-1))^2;
                end
            case 'EG2'
                for i = 1:n_obj-1
                    f = f+ sin(w(i)+w(i)^2-1)+0.5*sin(w(i)^2);
                end
            case 'Edensch'
                f = f+16;
                for i = 1:n_obj-1
                    f = f+ (w(i)-2)^4+(w(i)*w(i+1)-2*w(i+1))^2+(w(i+1)+1)^2;
                end
            case 'Indef'
                f  = f+sum(w);
                for i = 2:n_obj-1
                    f = f+ 0.5*cos(2*w(i)-w(end)-w(1));
                end
            case 'Cube'
                f = f+(w(1)-1)^2;
                for i = 2:n_obj
                    f = f+ 100*(w(i)-w(i-1)^3)^2;
                end
            case 'Bdexp'
                for i = n_obj-2
                    f = f+(w(i)-w(i+1))*exp(-w(i+2)*(w(i)+w(i+1)));
                end
            case 'Genhumps'
                for i = n_obj-1
                    f = f+sin(2*w(i))^2*sin(2*w(i+1))^2+0.05*(w(i)^2+w(i+1)^2);
                end
            case 'Arwhead'
                for i = n_obj-1
                    f = f + (-4*w(i)+3)+(w(i)^2+w(end)^2)^2;
                end
            case 'Quartc'
                for i = n_obj
                    f = f + (w(i)-1)^4;
                end
            case 'Cosine'
                for i = n_obj-1
                    f = f + cos(-0.5*w(i+1)+w(i)^2);
                end
            case 'Sine'
                for i = n_obj-1
                    f = f + sin(-0.5*w(i+1)+w(i)^2);
                end

                % last four are from the cutest https://ccom.ucsd.edu/~optimizers/cutest/problems/
            case 'LUKVLE10'
                f = 0;
                for i = 1:(n_obj/2)
                    xi = w(2*i-1);
                    yi = w(2*i);
                    f = f + (xi^2)^(yi^2 + 1) + (yi^2)^(xi^2 + 1);
                end
            case 'CURLY20'
                K = 20;
                f = 0;
                for i = 1:n_obj
                    if i <= n_obj-K
                        q = sum(w(i:(i+K)));
                    else
                        q = sum(w(i:n_obj));
                    end
                    f = f + q * (q * (q^2 - 20) - 0.1);
                end
            case 'SCURLY30'
                if n_obj > 30
                    K = 30;
                else
                    K = n_obj;
                end

                % SCAL = 12.0;
                SCAL = 0.01;
                f = 0;
                % scale factors as column vector
                S = exp(((0:(n_obj-1))/(n_obj-1)) * SCAL)';
                for i = 1:n_obj
                    if i <= n_obj-K
                        q = sum(S(i:(i+K)) .* w(i:(i+K)));
                        % else
                        %     q = sum(S(i:end) .* w(i:end));
                    end
                    f = f + q * (q * (q^2 - 20) - 0.1);
                end

            case 'SCURLY10'
                scal = 12.0;
                if n_obj > 10
                    k = 10;
                else
                    k = n_obj;
                end

                % Compute scaling factors S(i)
                S = zeros(n_obj, 1);
                for i = 1:n_obj
                    rat = (i-1) / (n_obj-1);
                    S(i) = exp(rat * scal);
                end

                % Compute objective function
                for i = 1:min(n_obj-k, n_obj)
                    % Sum from j=i to min(i+k, n_obj)
                    gvar = 0;
                    for j = i:min(i+k, n_obj)
                        gvar = gvar + S(j) * w(j);
                    end

                    % Apply P4 group function: gvar * (gvar * (gvar^2 - 20) - 0.1)
                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end

                % Handle remaining groups for i > n_obj-k
                for i = max(n_obj-k+1, 1):n_obj
                    gvar = 0;
                    for j = i:n_obj
                        gvar = gvar + S(j) * w(j);
                    end

                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end

            case 'SCURLY20'
                if n_obj > 20
                    k = 20;
                else
                    k = n_obj;
                end
                scal = 12.0;

                % Compute scaling factors S(i)
                S = zeros(n_obj, 1);
                for i = 1:n_obj
                    rat = (i-1) / (n_obj-1);
                    S(i) = exp(rat * scal);
                end

                % Compute objective function
                for i = 1:min(n_obj-k, n_obj)
                    % Sum from j=i to min(i+k, n_obj)
                    gvar = 0;
                    for j = i:min(i+k, n_obj)
                        gvar = gvar + S(j) * w(j);
                    end

                    % Apply P4 group function: gvar * (gvar * (gvar^2 - 20) - 0.1)
                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end

                % Handle remaining groups for i > n_obj-k
                for i = max(n_obj-k+1, 1):n_obj
                    gvar = 0;
                    for j = i:n_obj
                        gvar = gvar + S(j) * w(j);
                    end

                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end


            case 'CURLY10'
                 if n_obj > 10
                    k = 10;
                else
                    k = n_obj;
                end

                % Compute objective function
                for i = 1:min(n_obj-k, n_obj)
                    % Sum from j=i to min(i+k, n_obj)
                    gvar = 0;
                    for j = i:min(i+k, n_obj)
                        gvar = gvar + w(j);
                    end

                    % Apply P4 group function: gvar * (gvar * (gvar^2 - 20) - 0.1)
                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end

                % Handle remaining groups for i > n_obj-k
                for i = max(n_obj-k+1, 1):n_obj
                    gvar = 0;
                    for j = i:n_obj
                        gvar = gvar + w(j);
                    end

                    f = f + gvar * (gvar * (gvar^2 - 20) - 0.1);
                end
            case 'FMINSURF'
                p = round(sqrt(n_obj));  % assuming n_obj = p^2
                p_minus_1 = max(p - 1,1);
                scale = 1 / (p_minus_1)^2;
                param = 0.5 * (p_minus_1)^2;

                % Surface area terms (SQROOT groups)
                for i = 1:p_minus_1
                    for j = 1:p_minus_1
                        % Convert 2D indices to 1D
                        idx_ij = (j-1)*p + i;
                        idx_i1j1 = j*p + (i+1);
                        idx_i1j = (j-1)*p + (i+1);
                        idx_ij1 = j*p + i;

                        % a(i,j) = x(i,j) - x(i+1,j+1)
                        a_ij = w(idx_ij) - w(idx_i1j1);
                        % b(i,j) = x(i+1,j) - x(i,j+1)
                        b_ij = w(idx_i1j) - w(idx_ij1);

                        alpha = 1 + param * (a_ij^2 + b_ij^2);
                        f = f + scale * sqrt(alpha);
                    end
                end

                % Average height term (L2 group)
                avg_height = 0;
                for idx = 1:n_obj
                    avg_height = avg_height + w(idx);
                end
                f = f + (p^4) * scale * avg_height^2;

            case 'FLETBV3M'
                f = 0;
                for i = 1:n_obj
                    % finite-difference like group (difference with next)
                    if i < n_obj
                        q = w(i+1) - w(i);
                    else
                        q = w(i);
                    end
                    % element contribution (scaled sin/cos)
                    f = f + 0.5 * q^2 + 100*sin(0.01*w(i)) + 1e2*cos(w(i));
                end

            case 'SPARSQUR'
                for i = 1:n_obj
                    % Calculate indices based on SIF modular arithmetic
                    j2 = mod(2*i - 1, n_obj) + 1;
                    j3 = mod(3*i - 1, n_obj) + 1;
                    j5 = mod(5*i - 1, n_obj) + 1;
                    j7 = mod(7*i - 1, n_obj) + 1;
                    j11 = mod(11*i - 1, n_obj) + 1;

                    % Sum of squares with weight i
                    alpha = w(i) + w(j2) + w(j3) + w(j5) + w(j7) + w(j11);
                    f = f + 0.5 * i * alpha^2;
                end

            case 'NCVXQP6'

                % number of positive/negative rank-one terms
                nplus = floor(0.75 * n_obj);

                % Objective: quadratic form with sign-flipped rank-one terms
                f = 0;
                for i = 1:n_obj
                    % each group uses w(i), w(mod(2i-1,n)+1), w(mod(3i-1,n)+1)
                    j = mod(2*i-1, n_obj) + 1;
                    k = mod(3*i-1, n_obj) + 1;

                    alpha = w(i) + w(j) + w(k);

                    if i <= nplus
                        f = f + 0.5 * (alpha^2);    % positive term
                    else
                        f = f - 0.5 * (alpha^2);    % negative term
                    end
                end
            case 'DIXCHLNV'
                % w = decision vector of length n_obj
                % SIF requires n_obj even and >= 4
                f = 0;
                for i = 1:(n_obj-3)
                    % elements (from SIF):
                    % XSQ(k) = w(k)^2 for k = 1..n-1
                    xsq_i   = w(i)^2;           % XSQ(i)
                    xsq_ip1 = w(i+1)^2;         % XSQ(i+1)
                    xsq_ip2 = w(i+2)^2;         % XSQ(i+2)
                    xsq_ip3 = w(i+3)^2;         % XSQ(i+3)
                    % PR(i) = (w(i+1)-1)*(w(i+3)-1)
                    pr_i = (w(i+1)-1)*(w(i+3)-1);

                    % group A: scale 0.01, uses -XSQ(i)  -> (0.01 * (-xsq_i))^2
                    f = f + (0.01 * (-xsq_i))^2;

                    % group B: (unscaled) use xsq_i  -> (xsq_i)^2
                    f = f + (xsq_i)^2;

                    % group C: scale 1/90 on x(i+3) element -> ((1/90)*xsq_ip3)^2
                    f = f + ((1/90.0) * xsq_ip3)^2;

                    % group D: use xsq(i+2)
                    f = f + (xsq_ip2)^2;

                    % group E: scale 1/10.1 on x(i+1) -> ((1/10.1)*xsq_ip1)^2
                    f = f + ((1/10.1) * xsq_ip1)^2;

                    % group F: scale 1/10.1 on x(i+3) -> ((1/10.1)*xsq_ip3)^2
                    f = f + ((1/10.1) * xsq_ip3)^2;

                    % group G: scale 1/19.8 on PR(i) -> ((1/19.8)*pr_i)^2
                    f = f + ((1/19.8) * pr_i)^2;
                end
            case 'EDENSCH'
                for i = 1:n_obj-1
                    r = w(i)*w(i+1);        % pairwise product
                    f = f + r^2 + r^4;      % L2 + L4 contributions
                end


                % case 'Hager'
                %very similar to  Raydan 2
                %     for i = 1:n_obj-1
                %         f = f+(exp(w(i))-sqrt(i)*w(i));
                %     end
                % case 'Diagonal1'
                %  %very similar to  Raydan 2
                %     for i = 1:n_obj
                %         f = f+(exp(w(i))-i*w(i));
                %     end
                % case 'Diagonal2'
                % very similar to  Raydan 1
                %     for i = 1:n_obj-1
                %         f = f+(exp(w(i))-w(i)/i);
                %     end
        end

        f = f*1e-3; % alliviate bad scaling

        % Generate problem matrices
        if settings.inequality_constraints
            % Ax + By >= f
            s_density_A = s_density_A_B;
            if settings.nnz_bounded_by_dim 
                if ~(n_comp <= n_comp_min && adaptive_density_bounds)
                    s_density_A = nnz_factor*n_ineq/(n_x*n_ineq);
                end
            end
            if isempty(inv_cond_num)
                A = (range_A(2)-range_A(1))*sprand(n_ineq,n_x,s_density_A);
            else
                A = (range_A(2)-range_A(1))*sprand(n_ineq,n_x,s_density_A,inv_cond_num);
            end
    
            A(A~=0) = range_A(1)+A(A~=0);
            A = round(A,n_digits_data);
            if settings.inequality_constraints_coupling_terms
                s_density_B = s_density_A_B;
                if settings.nnz_bounded_by_dim
                    if ~(n_comp <= n_comp_min && adaptive_density_bounds)
                        s_density_B = nnz_factor*n_ineq/(n_y*n_ineq);
                    end
                end
                if isempty(inv_cond_num)
                    B = (range_B(2)-range_B(1))*sprand(n_ineq,n_y,s_density_B);
                else
                    B = (range_B(2)-range_B(1))*sprand(n_ineq,n_y,s_density_B ,inv_cond_num);
                end

                B(B~=0) = range_B(1)+B(B~=0);
                B = round(B,n_digits_data);
            else
                B = zeros(n_ineq,n_y);
            end
        end

        r = round(1+(n_y-1)*rand(1)); % pick number between
        if r == n_y
            r = r-1;
        end

         
         
        % s_density_M = (n_non_zero_E-n_y)/n_y^2;
        if settings.variable_density
            % s_density_M = range_s_density(1)+(range_s_density(2)-range_s_density(1)).*rand(1);
            s_density_M = range_s_density_lb +(range_s_density_ub-range_s_density_lb).*rand(1);
            s_density_M = min(s_density_M, n_non_zero_E/(r*(n_y-r))); % given desnity but not more than 2000 nnz
            if settings.nnz_bounded_by_dim
                if ~(n_comp <= n_comp_min && adaptive_density_bounds)
                    s_density_M = nnz_factor*min(r,n_non_zero_E)/(n_y*n_ineq);
                end
            end
        else
            s_density_M = settings.s_density_M;
        
        end
        if isempty(inv_cond_num)
            E = (range_E(2)-range_E(1))*sprand(r,n_y-r,s_density_M);
        else
            E = (range_E(2)-range_E(1))*sprand(r,n_y-r,s_density_M,inv_cond_num);
        end

        E(E~=0) = range_E(1)+E(E~=0);

        d1 = range_d1(1)+(range_d1(2)-range_d1(1)).*rand(r,1);
        d2 = range_d2(1)+(range_d2(2)-range_d2(1)).*rand(n_y-r,1);
        D1 = diag(d1);
        D2 = diag(d2);
        if settings.symmetric_psd_M_matrix
            M = [D1 E;...
                -E' D2];
            [V,D] = eig(full(M));
            M = V*abs(D)*V';
        else
            M = [D1 E;...
                -E' D2];
        end
        % M = (range_A(2)+1)*sprand(n_y,n_y,s_density_A_B);
        % M (M ~=0) = -1+M(M~=0);

        M = round(M,n_digits_data);

        % N = sprand(n_y,n_x,s_density);
        % r = a + (b-a).*rand(100,1);
        s_density_N = s_density_A_B;
        if settings.nnz_bounded_by_dim
            if ~(n_comp <= n_comp_min && adaptive_density_bounds)
                s_density_N = nnz_factor*n_y/(n_x*n_y);
            end
        end
        if isempty(inv_cond_num)
            N = (range_N(2)-range_N(1))*sprand(n_y,n_x,s_density_N);
        else
            N = (range_N(2)-range_N(1))*sprand(n_y,n_x,s_density_N,inv_cond_num);
        end

        N(N~=0) = range_N(1)+N(N~=0);
        N = round(N,n_digits_data);

        q = range_q(1)+(range_q(2)-range_q(1))*rand(n_y,1);
        q = round(q,n_digits_data);

        % feasible_ineq;
        if settings.inequality_constraints
            if settings.use_normal_distributions_for_feasible_point
                x_bar = abs(normrnd(0,1,[n_x 1]));
                % x_bar(x_bar<1) = 0;
                y_bar = max(0,normrnd(0,1,[n_y 1]));
                % y_bar = max(0,normrnd(-0.5,1,[n_y 1]));
            else
                x_bar = range_x(1)+(range_x(2)-range_x(1)).*rand(n_x,1);
                y_bar = range_y(1)+(range_y(2)-range_y(1)).*rand(n_y,1);
            end
            f_ineq_eps = abs(normrnd(0,1,[n_ineq, 1]));
            % f_ineq_eps  = 0;
            f_ineq = A*x_bar+B*y_bar-f_ineq_eps;
            f_ineq = round(f_ineq,n_digits_data);
        end
        % Inequality constraints
        if settings.inequality_constraints
            g_ineq = A*x+B*y-f_ineq;
            lbg_ineq = zeros(n_ineq,1);
            ubg_ineq = inf*ones(n_ineq,1);
        else
            g_ineq = [];
            lbg_ineq = [];
            ubg_ineq = [];
        end

        if nonlinear_ineq_constraints && ~strcmp(objective_type,'Linear') && ~strcmp(objective_type,'Quadratic_psd') && ~strcmp(objective_type,'Quadratic_ind')
            ind_nonlinear_ineq  = sort(randperm(n_ineq,round(settings.s_nonlinear_ineq*n_ineq)));
            if ~isempty(ind_nonlinear_ineq)
                g_ineq(ind_nonlinear_ineq) =  g_ineq(ind_nonlinear_ineq)+g_ineq(ind_nonlinear_ineq).^2+g_ineq(ind_nonlinear_ineq).^4;
            end
        end
        %
        if settings.copy_ineq
            ind_ineq_copy  = sort(randperm(n_ineq,round(settings.s_ineq_copy*n_ineq)));
            g_ineq = [g_ineq;g_ineq(ind_ineq_copy)];
            lbg_ineq = [lbg_ineq;lbg_ineq(ind_ineq_copy)];
            ubg_ineq = [ubg_ineq;ubg_ineq(ind_ineq_copy)];
        end

        % Equality constraints (lifted complementarity constraints)
        g_eq = N*x+M*y-z;
        if nonlinear_eq_constraints && ~strcmp(objective_type,'Linear') && ~strcmp(objective_type,'Quadratic_psd') && ~strcmp(objective_type,'Quadratic_ind')
            ind_nonlinear_eq  = sort(randperm(n_y,round(settings.s_nonlinear_eq*n_y)));
            if ~isempty(ind_nonlinear_eq)
                g_eq(ind_nonlinear_eq) =  g_eq(ind_nonlinear_eq)+g_eq(ind_nonlinear_eq).^2+g_eq(ind_nonlinear_eq).^4;
            end
        end
        lbg_eq = zeros(n_y,1);
        ubg_eq = zeros(n_y,1);

        % Copy equality constraints
        if settings.copy_eq
            ind_eq_copy  = sort(randperm(n_y,round(settings.s_eq_copy*n_y)));
            g_eq = [g_eq;g_eq(ind_eq_copy)];
            lbg_eq = [lbg_eq;lbg_eq(ind_eq_copy)];
            ubg_eq = [ubg_eq;ubg_eq(ind_eq_copy)];
        end



        n_eq = length(g_eq);
        n_ineq = length(g_ineq);

        lbw = [zeros(n_x+n_y,1);-inf*ones(n_y,1)];
        if bounded_w
            ubw = range_ubw(1)+(range_ubw(2)-range_ubw(1)).*rand(n_x+n_y,1);
            ubw = round(ubw,n_digits_data);
            ubw = [ubw;inf*ones(n_y,1)];
        else
            ubw = inf*ones(n,1);
        end
        % create problem functions
        % f = f;
        G = y;
        H = z;
        g = [g_eq;g_ineq];
        lbg = [lbg_eq;lbg_ineq];
        ubg = [ubg_eq;ubg_ineq];

        G_fun = Function('G_fun',{w},{G});
        H_fun = Function('H_fun',{w},{H});
        f_fun = Function('f_fun',{w},{f});
        g_fun = Function('g_fun',{w},{g});

        name = [problem_set_name '_' objective_type '_Var' num2str(n) '_Comp' num2str(n_comp) '_Eq' num2str(n_eq) '_Ineq' num2str(n_ineq)];
        w0 = [x_bar;y_bar;q+N*x_bar+M*y_bar];
        % save problem
        mpec.w = w;
        mpec.f_fun = f_fun;
        mpec.g_fun = g_fun;
        mpec.G_fun = G_fun;
        mpec.H_fun = H_fun;
        mpec.f = f;
        mpec.g = g;
        mpec.G = G;
        mpec.H = H;
        mpec.w0 = w0;
        mpec.lbw = lbw;
        mpec.ubw = ubw;
        mpec.lbg = lbg;
        mpec.ubg = ubg;
        mpec.name = name;
        mpecs = [mpecs, mpec];
        % end
    end
end
% save('NonLinMPECdat','mpecs')
end

