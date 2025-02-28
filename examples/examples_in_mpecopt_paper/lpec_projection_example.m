clear; clc; close all

%% An attemp to make the projection lpee fail
import casadi.*
a = 100;
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
f_sym = @(x,y) -2*x+y;
a = 0.1;
a = 1e-2;
g_sym  = @(x,y) y-(x)^2-a; % good lpec inconsistent bnlp
g_nabla_sym = @(x,y) [-2*x;1];
% g_sym  = @(x,y) y-(x)^2; % hard to verify B
f = f_sym(x1,x2);
g = [g_sym(x1,x2)];
x0 = [0.0;0.5];

lbx = [0;0];
ubx = [inf;1];
G = x1;
H = x2;
lbg = 0.0;
ubg = inf;

mpec = struct('x', x,'f',f, 'g',g,'G',G,'H',H);
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
% homotopy
settings = HomotopySolverOptions();
settings.tol  = 1e-12;
settings.max_iter = 15;
[result_homotopy,stats_homotopy] = mpec_homotopy_solver(mpec,solver_initalization,settings);
f_opt_homotopy = full(result_homotopy.f);
w_opt_homotopy = full(result_homotopy.x);
%%  Settings
% x0 = [1/a;0];
% x0 = [2;0];
solver_settings = mpecopt.Options();
solver_settings.settings_lpec.lpec_solver ="Gurobi";
solver_settings.initialization_strategy = "RelaxAndProject";
x0 = [0;0.5];
solver_initalization = struct('x0', x0, 'lbx',lbx, 'ubx',ubx,'lbg',lbg, 'ubg',ubg);
[result_active_set,stats_active_set] = mpec_optimizer(mpec, solver_initalization, solver_settings);
w_opt_active_set = full(result_active_set.x);
f_opt_active_set = full(result_active_set.f);

fprintf('x_homotopy = (%2.4f,%2.4f), f_opt = %2.4f. \n',w_opt_homotopy(1),w_opt_homotopy(2),f_opt_homotopy);
fprintf('x_active_set = (%2.4f,%2.4f), f_opt = %2.4f.  \n',w_opt_active_set(1),w_opt_active_set(2),f_opt_active_set);

%%
nice_plot_colors
tt = 0:1:5;
figure
plot(tt,tt*0,'k','LineWidth',1.5);
hold on
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
axis equal
plot(stats_active_set.iter.X_outer(1,:),stats_active_set.iter.X_outer(2,:),'rs')
plot(stats_active_set.iter.X_lpec(1,:),stats_active_set.iter.X_lpec(2,:),'bd')
if 1
    for i = 1:size(stats_active_set.iter.X_outer,2)-1
        quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), stats_active_set.iter.X_outer(1,i+1)-stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i+1)-stats_active_set.iter.X_outer(2,i), 1, 'Color',matlab_red,'LineWidth',2);
        text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
    end
else
    for i = 1:size(stats_active_set.iter.X_outer,2)-1
        quiver(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i), -2*(a*stats_active_set.iter.X_outer(1,i)-1),-2*(stats_active_set.iter.X_outer(2,i)-1), 0.5, 'Color',matlab_red,'LineWidth',2);
        text(stats_active_set.iter.X_outer(1,i), stats_active_set.iter.X_outer(2,i),num2str(i));
    end
end
xlim([-1 5])
ylim([-1 5])

f = @(x,y) (a*x-1).^2+(y-1).^2;
% Generate grid points for x and y
x = linspace(-1, 5, 50);
y = linspace(-1, 5, 50);
% Create a grid of points
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f_sym(X, Y);
% Plot the contour lines
contour(X, Y, Z, 30);
fimplicit(g_sym, [-2 2])
xlabel('x');
ylabel('y');
%%
filename1  = 'cq_mpec_1.pdf';
filename2  = 'cq_lpec_2.pdf';
rho = 0.5;
nice_plot_colors;
linewidht = 2;
fontsize = 16;
markersize = 10;
x_opt = [0;0];
nice_plot_colors
tt = 0:1:5;
figure
latexify_plot();
x = linspace(-2, 5, 50);
y = linspace(-2, 5, 50);
[X, Y] = meshgrid(x, y);
% Evaluate the function at each point in the grid
Z = f_sym(X, Y);
% Plot the contour lines
contour(X, Y, Z, 30);
hold on;
fimplicit(g_sym, [-5 5], 'MeshDensity', 100, 'LineWidth',2,'Color',matlab_blue)

plot(tt,tt*0,'k','LineWidth',1.5);
hold on
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
axis equal
plot(x_opt(1),x_opt(2),'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
text(x_opt(1)-0.4,x_opt(2)+0.3, '${x}^*$','FontSize',fontsize)

xlim([-1.2 3])
ylim([-1.2 3])
xlabel('$x_1$');
ylabel('$x_2$');

set(gca,'FontSize',fontsize);
legend({'$f(x)$', '$c(x)=0$' },'BackgroundAlpha',0.9)
exportgraphics(gcf, filename1, 'ContentType', 'vector')

% Linearization
figure
latexify_plot();
contour(X, Y, Z, 30);
hold on
quiver(x_opt(1),x_opt(2),2,-1,'LineWidth',linewidht+0.5,'AutoScaleFactor',1,'Color',matlab_red)
fimplicit(g_sym, [-5 5], 'MeshDensity', 100, 'LineWidth',1, 'LineStyle','--','Color',matlab_blue)
fimplicit(@(x,y) g_sym(x_opt(1),x_opt(2))+g_nabla_sym(x_opt(1),x_opt(2))'*([x;y]-x_opt), [-5 5], 'MeshDensity', 100, 'LineWidth',linewidht,'Color',matlab_blue)
% optimal lpec solution
% quiver(x_opt(1),x_opt(2),-x_opt(1)+0.5,-x_opt(2)-0.1,'LineWidth',linewidht+0.5,'Color',matlab_magenta,'MaxHeadSize',0.8)
quiver(x_opt(1),x_opt(2),rho,0,'LineWidth',linewidht+0.5,'Color',matlab_magenta,'MaxHeadSize',0.8)


% Trust region radius
x = x_opt;
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
hold on
plot(tt,tt*0-rho+x(2),'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',matlab_blood_red,'LineWidth',linewidht+0.5,'HandleVisibility','off');


tt = 0:1:5;
plot(tt,tt*0,'k','LineWidth',1.5);
grid on
plot(tt*0,tt,'k','LineWidth',1.5);
% overwrite 
fimplicit(@(x,y) g_sym(x_opt(1),x_opt(2))+g_nabla_sym(x_opt(1),x_opt(2))'*([x;y]-x_opt), [-5 5], 'MeshDensity', 100, 'LineWidth',linewidht,'Color',matlab_blue)
quiver(x_opt(1),x_opt(2),rho,0,'LineWidth',linewidht+0.5,'Color',matlab_magenta,'MaxHeadSize',0.8)
plot(x_opt(1),x_opt(2),'Marker','pentagram','MarkerEdgeColor',matlab_orange,'MarkerFaceColor',matlab_orange,'MarkerSize',markersize)
text(x_opt(1)-0.4,x_opt(2)+0.3, '${x}^*$','FontSize',fontsize)

axis equal
xlim([-1.2 3])
ylim([-1.2 3])
legend({'$f(x)$', '$-\nabla f(x^*(\tau))$', '$c(x)=0$', '$c(x^*(\tau))\!+\! \nabla c(x^*(\tau))^\top d=0$', '$d$' },'BackgroundAlpha',0.8,'NumColumns',1,'FontSize',fontsize-2)
xlabel('$x_1$');
ylabel('$x_2$');
set(gca,'FontSize',fontsize);

exportgraphics(gcf, filename2, 'ContentType', 'vector')

