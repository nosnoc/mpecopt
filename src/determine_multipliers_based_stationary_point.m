function [multiplier_based_stationarity, output_message] = determine_multipliers_based_stationary_point(x,lambda,dims,settings)


eta_x1  = -lambda(dims.ind_x1); % MPEC Lagrange multiplieres of x1\geq (=) 0
eta_x2  = -lambda(dims.ind_x2); % MPEC Lagrange multiplieres of x2 \geq (=) 0


eta_x1(abs(eta_x1) < settings.tol*10) = 0;
eta_x2(abs(eta_x2) < settings.tol*10) = 0;
eta_x1(abs(eta_x1) < settings.tol) = 0;
eta_x2(abs(eta_x2) < settings.tol) = 0;

active_set_estimate_k = find_active_sets(x, dims, settings.tol_active);
ind_I_00 = find(active_set_estimate_k.I_00);

% Isolate multipliers for biactives

eta_x1_biactive = eta_x1(ind_I_00);
eta_x2_biactive = eta_x2(ind_I_00);
max_mult = max(norm([eta_x1_biactive;eta_x2_biactive],inf), 1e-12);
% max_mult = 1;

% Rescale to max (sensitive for M-stationarity)
eta_x1_biactive = eta_x1_biactive./max_mult; 
eta_x2_biactive = eta_x2_biactive./max_mult; 

% if max_mult > 1e3
%     tol_M = 1e-3;
% else
tol_M = settings.tol_B_stationarity*10;
% end


% R

% print_iter_line();
if isempty(ind_I_00)
    output_message = sprintf('S - stationary, %d biactive constraints. ', length(ind_I_00));
    multiplier_based_stationarity = 'S';
elseif (all(eta_x1_biactive >= 0 & eta_x2_biactive >= 0) )
    output_message = sprintf('S - stationary, %d biactive constraints. ', length(ind_I_00));
    multiplier_based_stationarity = 'S';
    if all(eta_x1.*eta_x2 > 0)
        output_message = [output_message ,' Upper level strict complementarity holds.'];
    else
        output_message = [output_message ,' Upper level strict complementarity does not hold.'];
    end
else
    % if all( (eta_x1_biactive > 0 & eta_x2_biactive > 0) | (eta_x1_biactive.*eta_x2_biactive == 0))
    if all( (eta_x1_biactive > 0 & eta_x2_biactive > 0) | abs(eta_x1_biactive.*eta_x2_biactive)<= tol_M)
        % account for imperfect zeros
        output_message = sprintf('M - stationary, %d biactive constraints.', length(ind_I_00));
        multiplier_based_stationarity = 'M';
    elseif all(eta_x1_biactive.*eta_x2_biactive >=0)
        output_message = sprintf('C - stationary, %d biactive constraints.', length(ind_I_00));
        multiplier_based_stationarity = 'C';
    elseif  all(eta_x1_biactive >= 0 | eta_x2_biactive >=0)
        output_message = sprintf('A - stationary, %d biactive constraints.', length(ind_I_00));
        multiplier_based_stationarity = 'A';
    else
        output_message = sprintf('W - stationary, %d biactive constraints.', length(ind_I_00));
        multiplier_based_stationarity = 'W';
    end
end
if settings.plot_mpec_multipliers
    if ~isempty(ind_I_00)
        % figure
        % clf;
        hold on;
        switch multiplier_based_stationarity
            case 'S'
                title('S')
                hh_x = [0  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 0];
                hh_y = [0 0 max(abs(eta_x2(:)))+1 max(abs(eta_x2(:)))+1 ];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
                % keyboard
            case 'M'
                title('M')
                hh_x = [0  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 0];
                hh_y = [0 0 max(abs(eta_x2(:)))+1 max(abs(eta_x2(:)))+1 ];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
                xline(0,'LineWidth',1.5,'Color',[0 0 0 0.1])
                yline(0,'LineWidth',1.5,'Color',[0 0 0 0.1])
            case 'C'
                title('C')
                hh_x = [0  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 0];
                hh_y = [0 0 max(abs(eta_x2(:)))+1 max(abs(eta_x2(:)))+1 ];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
                hh_x = [0  -max(abs(eta_x1(:)))-1  -max(abs(eta_x1(:)))-1 0];
                hh_y = [0 0 -max(abs(eta_x2(:)))-1 -max(abs(eta_x2(:)))-1 ];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
            case 'A'
                title('A')
                hh_x = [-max(abs(eta_x1(:)))-1  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 -max(abs(eta_x1(:)))-1];
                hh_y = [0  0  max(abs(eta_x2(:)))+1 max(abs(eta_x2(:)))+1];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;

                hh_x = [0  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 0];
                hh_y = [0 0 -max(abs(eta_x2(:)))-1 -max(abs(eta_x2(:)))-1];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
            case 'W'
                title('W')
                hh_x = [-max(abs(eta_x1(:)))-1  max(abs(eta_x1(:)))+1  max(abs(eta_x1(:)))+1 -max(abs(eta_x1(:)))-1];
                hh_y = [-max(abs(eta_x2(:)))-1  -max(abs(eta_x2(:)))-1  max(abs(eta_x2(:)))+1 max(abs(eta_x2(:)))+1];
                hh = patch(hh_x,hh_y,'k');
                hh.FaceAlpha = 0.1;
        end
        plot(eta_x1_biactive,eta_x2_biactive,'.','LineWidth',2.5,'MarkerSize',15)
        grid on
        xlabel('$\nu$','Interpreter','latex')
        ylabel('$\xi$','Interpreter','latex')
        hold on
        xline(0,'k')
        hold on
        yline(0,'k')
        xlim([-max(abs(eta_x1(:)))-1 max(abs(eta_x1(:)))+1])
        ylim([-max(abs(eta_x2(:)))-1 max(abs(eta_x2(:)))+1])
        %                             axis equal
    end
end

end