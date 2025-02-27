clear all; clc; close all;
%%
import casadi.*
% different stationarity types
origin_stationarity = 'M';
switch origin_stationarity 
    case 'S'
%         a = 1;  c = 0;  b = -1;  d =0 ;
        a = 1; c = 1; b = -1; d = -1;
    case 'M'    
        a = 1;  c = 0;  b = +1;  d =0 ;
%           a = 0;  c = 1;  b = 0;  d = 1 ;
    case 'C'
         a = 1;  c = 1; b = +1; d =+1;
    case 'A'
        a = 1;  c = 1;  b = -1;  d = +1 ;
        % a = 0;  c = 1;  b = 0;  d = 1 ;
end

x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];

f = (a*x1-b)^2+(c*x2-d)^2;
G = x1;
H = x2;

% x0 = [2;0];
x0  = [0;0];

lbx = [0;0];
ubx = [inf;inf];
g = []; lbg = []; ubg = [];
mpec = struct('x', x, 'f', f, 'g', g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
%% reg
settings = HomotopySolverOptions();
[sol,stats] =  mpec_homotopy_solver(mpec,solver_initalization,settings);
w_opt_reg = full(sol.x);
f_opt_reg = full(sol.f);
%%  mpecopt
solver_settings = MPECOptimizerOptions();
solver_settings.initialization_strategy = "TakeInitialGuessDirectly";
% [result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);

solver = Mpecopt(mpec, solver_settings);
[result_active_set,stats_active_set] = solver.solve(solver_initalization);

w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);

fprintf('solution mpecopt is (%2.4f,%2.4f), f = %2.4f \n',w_opt_active_set(1),w_opt_active_set(2),f_opt_active_set);
fprintf('solution reg is (%2.4f,%2.4f), f = %2.4f \n',w_opt_reg(1),w_opt_reg(2),f_opt_reg );

%% Plot iterates
% nice_plot_colors
% figure
% xline(0,'k');
% hold on
% grid on
% yline(0,'k');
% axis equal
% plot(stats_active_set.iter.X_outer(1,:),stats_active_set.iter.X_outer(2,:),'rs')
% for i = 1:size(stats_active_set.iter.X_outer,2)-1
%     quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
%     text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
% end
% f = @(x1,x2) (a*x1-b).^2+(c*x2-d).^2;
% % Generate grid points for x and y
% quiver(0,0,-2*(a*0-b),-2*(c*0-d),'k','LineWidth',3);
% x = linspace(-1, 5, 50);
% y = linspace(-1, 5, 50);
% % Create a grid of points
% [X, Y] = meshgrid(x, y);
% % Evaluate the function at each point in the grid
% Z = f(X, Y);
% % Plot the contour lines
% contour(X, Y, Z, 10);
% xlabel('x');
% ylabel('y');

