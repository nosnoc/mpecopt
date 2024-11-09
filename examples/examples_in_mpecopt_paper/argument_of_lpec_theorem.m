%% Illustration of the argument in Theorem 4.1.
clear; close all;
nice_plot_colors
makersize = 10;
linewidht = 2.5;
fontsize = 16;
latexify_plot();
% plot at branch

x_opt = [1.0,0];
rho = 0.9;
x = x_opt;
tt = 0:0.1:3;
figure
plot(tt,tt*0,'k','LineWidth',linewidht,'HandleVisibility','off')
hold on
plot(tt*0,tt,'k','LineWidth',linewidht,'HandleVisibility','off')

% trust region radius  \B(x^*,\bar{\rho})
tr_color = matlab_green;
tt = linspace(x(1)-rho-1e-2,x(1)+rho+1e-2,3);
plot(x_opt(1),x_opt(2), 'Color',tr_color,'Marker','pentagram','MarkerFaceColor', tr_color,'MarkerEdgeColor',tr_color,'MarkerSize',makersize,'HandleVisibility','off');
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{B}(x^*,\bar{\rho})$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
text(x(1)+0.1,x(2)+0.2,'$x^*$','FontSize',fontsize,'Color',tr_color);


% the region N
tr_color = matlab_blood_red;
rho = 1.1;
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{N}^*=\mathcal{B}(x^*,\kappa)\cap X$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');

% trust region radius  \B(x^*,\bar{\rho})
tr_color = matlab_magenta;
x = [1.0,0.3];
rho = 0.5;
tt = linspace(x(1)-rho-1e-2,x(1)+rho+1e-2,3);
plot(x(1),x(2), 'Color',tr_color,'Marker','diamond','MarkerFaceColor', tr_color,'MarkerEdgeColor',tr_color,'MarkerSize',makersize-3,'HandleVisibility','off');
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{B}(\tilde{x},{\rho})$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
text(x(1)+0.1,x(2)+0.2,'$\tilde{x}$','FontSize',fontsize,'Color',tr_color);
% grid on
xlabel('$x_{1,i}$')
ylabel('$x_{2,i}$')
xticks([]);
xticklabels({});
yticks([]);
yticklabels({});


set(gca,"Fontsize",fontsize)
L = legend('BackgroundAlpha',0.7);
axis equal;
xlim([-0.5 3.0])
ylim([-1.5 2.0])
exportgraphics(gcf, 'proof_argument_1.pdf', 'ContentType', 'vector')


%% plot at corner
x_opt = [0.0,0];
rho = 1;
x = x_opt;
tt = 0:0.1:3;
figure
plot(tt,tt*0,'k','LineWidth',linewidht,'HandleVisibility','off')
hold on
plot(tt*0,tt,'k','LineWidth',linewidht,'HandleVisibility','off')

% trust region radius  \B(x^*,\bar{\rho})
tr_color = matlab_green;
tt = linspace(x(1)-rho-1e-2,x(1)+rho+1e-2,3);
plot(x_opt(1),x_opt(2), 'Color',tr_color,'Marker','pentagram','MarkerFaceColor', tr_color,'MarkerEdgeColor',tr_color,'MarkerSize',makersize,'HandleVisibility','off');
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{B}(x^*,\bar{\rho})$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
text(x(1)+0.1,x(2)+0.2,'$x^*$','FontSize',fontsize,'Color',tr_color);


% the region N
tr_color = matlab_blood_red;
rho = 1.2;
tt = linspace(x(1)-rho,x(1)+rho,3);
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{N}^*=\mathcal{B}(x^*,\kappa)\cap X$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');

% trust region radius  \B(x^*,\bar{\rho})
tr_color = matlab_magenta;
x = [0.4,0.3];
rho = 0.5;
tt = linspace(x(1)-rho-1e-2,x(1)+rho+1e-2,3);
plot(x(1),x(2), 'Color',tr_color,'Marker','diamond','MarkerFaceColor', tr_color,'MarkerEdgeColor',tr_color,'MarkerSize',makersize-3,'HandleVisibility','off');
plot(tt,tt*0+rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','on','DisplayName',"$\mathcal{B}(\tilde{x},{\rho})$");
hold on
plot(tt,tt*0-rho+x(2),'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
tt = linspace(x(2)-rho,x(2)+rho,3);
plot(tt*0+x(1)+rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
plot(tt*0+x(1)-rho,tt,'color',tr_color,'LineWidth',linewidht+0.5,'HandleVisibility','off');
text(x(1)+0.1,x(2)+0.2,'$\tilde{x}$','FontSize',fontsize,'Color',tr_color);
% grid on
xlabel('$x_{1,j}$')
ylabel('$x_{2,j}$')
xticks([]);
xticklabels({});
yticks([]);
yticklabels({});


set(gca,"Fontsize",fontsize)
L = legend('BackgroundAlpha',0.7);
axis equal;
xlim([-1.4 2.6])
ylim([-1.5 2.5])
exportgraphics(gcf, 'proof_argument_2.pdf', 'ContentType', 'vector')