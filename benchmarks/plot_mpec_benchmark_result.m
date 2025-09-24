function [] = plot_mpec_benchmark_result(dtable,filename,settings)

% latex stuff
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% add path for figures folder
filename = ['figures/' filename];

settings.use_reduced_lpec_data = 1; % make lpec plots only for lpecs, where at leas one was solved

%% todo plots;
% 1) stacked bar plot cpu time nlp per phase i and phase ii

% colors, style and order
linewidth = 2;
normal_figure_size = [100, 100, 600, 420];
wide_figure_size = [100, 100, 1100, 400];
color_order = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.6350 0.0780 0.1840
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330];

% color_order = repmat([linspace(0.2,1,7)]',1,3); % gray scale

line_style_cycling_method = 'withcolor';  % 'withcolor'; % 'beforecolor'; ' aftercolor'
% line_style_cycling_method = 'aftercolor';
line_style_order = ["-";"--";":";"-.";"--";":s"];  % ["-.";"-x";"-o";"-s"];
% line_style_order = ["-";"--";":";"-.";":x";":s"];  % ["-.";"-x";"-o";"-s"];
FontSize = 12;


% solver_tab = pivot(dtable, Rows= "solver_name");
% pivot does not keep original ordering, the next 3 lines fix this:
[~,temp_b] = unique(dtable.solver_name);
solver_name = dtable.solver_name(sort(temp_b));
solver_tab = table(solver_name);


dtable.cpu_time_lpec(dtable.cpu_time_lpec <1e-6) = 0; % fix that gurobi sets some times to zero.
dtable_ii = dtable(dtable.solver_name == solver_tab.solver_name(1),:);
ind_lpec = dtable_ii .n_lpec_total > 0;
dtable_lpec = [];
for ii = 1:height(solver_tab)
    dtable_ii = dtable(dtable.solver_name == solver_tab.solver_name(ii),:);
    dtable_lpec = [dtable_lpec;dtable_ii(ind_lpec,:)];
end
% solver_names

% dtable_lpecs


if settings.b_stationarty_as_success_criterion
    success_tab = dtable(dtable.success==true & dtable.b_stationarity==true, :);
    if settings.use_reduced_lpec_data
        success_tab_lpec = dtable_lpec(dtable_lpec.success==true & dtable_lpec.b_stationarity==true, :);
    else
        success_tab_lpec = success_tab;
    end
else
    success_tab = dtable(dtable.success==true, :);
    if settings.use_reduced_lpec_data
        success_tab_lpec = dtable_lpec(dtable_lpec.success==true , :);
    else
        success_tab_lpec = success_tab;
    end
end
success_tab_phase_i = dtable(dtable.success_phase_i==true, :);

minimum_times = pivot(success_tab, Rows= ["problem_name"], DataVariable="cpu_time", Method="min");
minimum_times_lpec = pivot(success_tab, Rows= ["problem_name"], DataVariable="cpu_time_lpec", Method="min");

minimum_times_phase_i = pivot(success_tab_phase_i, Rows= ["problem_name"], DataVariable="cpu_time_phase_i", Method="min");
minimum_times_phase_ii = pivot(success_tab, Rows = ["problem_name"], DataVariable="cpu_time_phase_ii", Method="min");

% calculate relative times
dtable.ratio_best = zeros(height(dtable),1);
dtable.ratio_best_phase_i = zeros(height(dtable),1);
dtable.ratio_best_phase_ii = zeros(height(dtable),1);
dtable_lpec.ratio_best_lpec = zeros(height(dtable_lpec),1);
for ii=1:height(dtable)
    try
        dtable.ratio_best(ii) = dtable.cpu_time(ii)/(minimum_times.min_cpu_time(minimum_times.problem_name == dtable.problem_name(ii)));
    catch
        dtable.ratio_best(ii) = nan;
    end

    try
        dtable.ratio_best_phase_i(ii) = dtable.cpu_time_phase_i(ii)/(minimum_times_phase_i.min_cpu_time_phase_i(minimum_times_phase_i.problem_name == dtable.problem_name(ii)));
    catch
        dtable.ratio_best_phase_i(ii) = nan;
    end

    try
        dtable.ratio_best_phase_ii(ii) = dtable.cpu_time_phase_ii(ii)/(minimum_times_phase_ii.min_cpu_time_phase_ii(minimum_times_phase_ii.problem_name == dtable.problem_name(ii)));
    catch
        dtable.ratio_best_phase_ii(ii) = nan;
    end

    try
        dtable_lpec.ratio_best_lpec(ii) = dtable_lpec.cpu_time_lpec(ii)/(minimum_times_lpec.min_cpu_time_lpec(minimum_times_lpec.problem_name == dtable_lpec.problem_name(ii)));
    catch
        dtable_lpec.ratio_best_lpec(ii) = nan;
    end
end

success_tab = dtable(dtable.success==true, :);
if settings.b_stationarty_as_success_criterion
    success_tab = dtable(dtable.success==true & dtable.b_stationarity==true, :);
    success_tab_lpec = dtable_lpec(dtable_lpec.success==true & dtable_lpec.b_stationarity==true, :);

else
    success_tab = dtable(dtable.success==true, :);
    success_tab_phase_i = dtable(dtable.success_phase_i==true, :);
    success_tab_lpec = dtable_lpec(dtable_lpec.success==true, :);
end



stationary_points_tab = pivot(dtable, Rows=["solver_name","multiplier_based_stationarity"]);
stationary_points_tab_with_b_stat = pivot(dtable, Rows=["solver_name","multiplier_based_stationarity","b_stationarity"]);

b_stationary_tab = pivot(dtable, Rows=["solver_name","b_stationarity"]);
% b_stationary_tab_among_success = pivot(success_tab, Rows=["solver_name","b_stationarity"]);

n_mpecs = length(dtable.cpu_time(dtable.solver_name == solver_tab.solver_name(1)));
n_mpecs_with_lpec_solves = length(dtable_lpec.cpu_time(dtable_lpec.solver_name == solver_tab.solver_name(1)));
step_size = 1/n_mpecs;
step_size_lpec = 1/n_mpecs_with_lpec_solves ;
%% Absolute plots
if settings.absolute
    handle_fig_absolute = figure('Position', normal_figure_size);
    xlabel('Wall time (s)')
    ylabel('Fraction solved')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_absolute);
        hold on;
        times = sort(success_tab.cpu_time(success_tab.solver_name == solver_name));
        levels = step_size*(1:length(times));
        label = [char(solver_name)];
        stairs([0,times',max(dtable.cpu_time)*1.1], [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off
    figure(handle_fig_absolute);
    set(gca,'xscale','log');
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    % legend('Location', 'southeast')
    legend('Location', 'east')
    % grid on;
    if settings.save_plot
        exportgraphics(gca, [filename '_absolute.pdf']);
    end
end

if settings.absolute_phase_i
    handle_fig_absolute = figure('Position', normal_figure_size);
    xlabel('Wall time (s)')
    ylabel('Fraction solved (Phase I)')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_absolute);
        hold on;
        times = sort(success_tab_phase_i.cpu_time_phase_i(success_tab_phase_i.solver_name == solver_name));
        levels = step_size*(1:length(times));
        label = [char(solver_name)];
        stairs([0,times',max(dtable.cpu_time_phase_i)*1.1], [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off
    figure(handle_fig_absolute);
    set(gca,'xscale','log');
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend('Location', 'southeast')
    
    % grid on;
    if settings.save_plot
        exportgraphics(gca, [filename '_absolute_phase_i.pdf']);
    end
end

if settings.absolute_phase_ii
    handle_fig_absolute = figure('Position', normal_figure_size);
    xlabel('Wall time (s)')
    ylabel('Fraction solved (Phase II)')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_absolute);
        hold on;
        times = sort(success_tab.cpu_time_phase_ii(success_tab.solver_name == solver_name));
        levels = step_size*(1:length(times));
        label = [char(solver_name)];
        stairs([0,times',max(dtable.cpu_time_phase_ii)*1.1], [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off
    figure(handle_fig_absolute);
    set(gca,'xscale','log');
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend('Location', 'southeast')
    % grid on;
    if settings.save_plot
        exportgraphics(gca, [filename '_absolute_phase_ii.pdf']);
    end
end


if settings.absolute_lpec
    handle_fig_absolute = figure('Position', normal_figure_size);
    xlabel('Total LPEC wall time (s)')
    ylabel('Fraction solved')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_absolute);
        hold on;
        times = sort(success_tab_lpec.cpu_time_lpec(success_tab_lpec.solver_name == solver_name));
        levels = step_size_lpec*(1:length(times));
        label = [char(solver_name)];
        stairs([0,times',max(dtable_lpec.cpu_time_lpec)*1.1], [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off
    figure(handle_fig_absolute);
    set(gca,'xscale','log');
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend('Location', 'southeast')
    % grid on;
    if settings.save_plot
        exportgraphics(gca, [filename '_absolute_lpec.pdf']);
    end
end

%% Relative plots
if settings.relative
    handle_fig_relative= figure('Position', normal_figure_size);
    xlabel('$2^{x}$ times best')
    ylabel('Fraction solved')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_relative);
        hold on;
        ratios = sort(success_tab.ratio_best(success_tab.solver_name == solver_name));
        if min(ratios) <1
            keyboard;
        end
        levels = step_size*(1:length(ratios));
        label = [char(solver_name)];
        if 1
            stairs(log2([0,ratios',2^15]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        else
            stairs(log10([0,ratios',10^5]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
        alpha(1-(ii-1)*0.1);
    end
    hold off
    figure(handle_fig_relative);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    % grid on;
    % legend('Location', 'southeast');
    legend('Location', 'east');
    if settings.save_plot
        exportgraphics(gca, [filename '_relative.pdf']);
    end
end


if settings.relative_phase_i
    handle_fig_relative= figure('Position', normal_figure_size);
    xlabel('$2^{x}$ times best')
    ylabel('Fraction solved (Phase I)')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_relative);
        hold on;
        ratios = sort(success_tab_phase_i.ratio_best_phase_i(success_tab_phase_i.solver_name == solver_name));
        if min(ratios) <1
            keyboard;
        end
        levels = step_size*(1:length(ratios));
        label = [char(solver_name)];
        if 1
            stairs(log2([0,ratios',2^15]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        else
            stairs(log10([0,ratios',10^5]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
        alpha(1-(ii-1)*0.1);
    end
    hold off
    figure(handle_fig_relative);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    % grid on;
    legend('Location', 'southeast');
    if settings.save_plot
        exportgraphics(gca, [filename '_relative_phase_i.pdf']);
    end
end


if settings.relative_phase_ii
    handle_fig_relative= figure('Position', normal_figure_size);
    xlabel('$2^{x}$ times best')
    ylabel('Fraction solved (Phase II)')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_relative);
        hold on;
        ratios = sort(success_tab.ratio_best_phase_ii(success_tab.solver_name == solver_name));
        levels = step_size*(1:length(ratios));
        label = [char(solver_name)];
        if 1
            stairs(log2([0,ratios',2^32]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        else
            stairs(log10([0,ratios',10^5]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
        alpha(1-(ii-1)*0.1);
    end
    hold off
    figure(handle_fig_relative);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    % grid on;
    legend('Location', 'southeast');
    if settings.save_plot
        exportgraphics(gca, [filename '_relative_phase_ii.pdf']);
    end
end


if settings.relative_lpec
    handle_fig_relative= figure('Position', normal_figure_size);
    xlabel('$2^{x}$ times best for LPEC')
    ylabel('Fraction solved')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_relative);
        hold on;
        ratios = sort(success_tab_lpec.ratio_best_lpec(success_tab_lpec.solver_name == solver_name));
        levels = step_size_lpec*(1:length(ratios));
        label = [char(solver_name)];
        if 1
            stairs(log2([0,ratios',2^32]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        else
            stairs(log10([0,ratios',10^5]), [0,levels,levels(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
        alpha(1-(ii-1)*0.1);
    end
    hold off
    figure(handle_fig_relative);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ylim([0, 1.0]);
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    % grid on;
    legend('Location', 'southeast');
    if settings.save_plot
        exportgraphics(gca, [filename '_relative_lpec.pdf']);
    end
end

%% Plot LPEC times per phase to see better effect of TR:
if settings.lpec_phases_cpu_time
    legend_str = {};
    all_data  = [];
    n_problems = length(dtable.cpu_time_lpec(dtable.solver_name==solver_tab.solver_name(1),:));
    cpu_max = 1.1*max(dtable.cpu_time_lpec_phase_i+dtable.cpu_time_lpec_phase_ii);
    for ii=1:height(solver_tab)
        handle_fig_bar_timings = figure('Position', normal_figure_size);
        % subplot(ceil(height(solver_tab)/2),2,ii)
        solver_name = solver_tab.solver_name(ii);
        hold on
        data_ii = [dtable.cpu_time_lpec_phase_i(dtable.solver_name==solver_name,:)';dtable.cpu_time_lpec_phase_ii(dtable.solver_name==solver_name,:)'];
        % data_ii = log10(data_ii);
        bar(1:n_problems,data_ii,'stacked')
        legend_str = [append(solver_name,' Phase I');append(solver_name,' Phase II')];
        legend(legend_str,'Location', 'northeast',BackgroundAlpha=0.5);
        xlabel('Problem instance')
        ylabel('Total LPEC wall time in seconds')
        ylim([0 cpu_max])
        figure(handle_fig_bar_timings);
        set(gca,'fontsize', FontSize);
        set(gca,'YScale','log');
        ax = gca;
        ax.ColorOrder= color_order;
        ax.LineStyleOrder = line_style_order;
        ax.LineStyleCyclingMethod = line_style_cycling_method;
        exportgraphics(gca, [filename '_bar_timing_' num2str(ii) '.pdf']);
    end
end


%% Plot NLP times per phase to see better effect of TR:
if settings.nlp_phases_cpu_time
    legend_str = {};
    all_data  = [];
    n_problems = length(dtable.cpu_time_nlp(dtable.solver_name==solver_tab.solver_name(1),:));
    cpu_max = 1.1*max(dtable.cpu_time_nlp_phase_i+dtable.cpu_time_nlp_phase_ii);
    for ii=1:height(solver_tab)
        handle_fig_bar_timings = figure('Position', normal_figure_size);
        % subplot(ceil(height(solver_tab)/2),2,ii)
        solver_name = solver_tab.solver_name(ii);
        hold on
        data_ii = [dtable.cpu_time_nlp_phase_i(dtable.solver_name==solver_name,:)';dtable.cpu_time_nlp_phase_ii(dtable.solver_name==solver_name,:)'];
        % data_ii = log10(data_ii);
        bar(1:n_problems,data_ii,'stacked')
        legend_str = [append(solver_name,' Phase I');append(solver_name,' Phase II')];
        legend(legend_str,'Location', 'northeast',BackgroundAlpha=0.5);
        xlabel('Problem instance')
        ylabel('Total NLP wall time in seconds')
        ylim([0 cpu_max])
        figure(handle_fig_bar_timings);
        set(gca,'fontsize', FontSize);
        set(gca,'YScale','log');
        ax = gca;
        ax.ColorOrder= color_order;
        ax.LineStyleOrder = line_style_order;
        ax.LineStyleCyclingMethod = line_style_cycling_method;
        exportgraphics(gca, [filename '_bar_timing_' num2str(ii) '.pdf']);
    end
end


%% Plot NLP and LPEC timings in bar plots
if settings.bar_timing_plots
    % handle_fig_bar_timings = figure('Position', normal_figure_size);
    legend_str = {};
    all_data  = [];
    n_problems = length(dtable.cpu_time_nlp(dtable.solver_name==solver_tab.solver_name(1),:));
    cpu_max = 1.1*max(dtable.cpu_time_nlp+dtable.cpu_time_lpec);
    for ii=1:height(solver_tab)
        handle_fig_bar_timings = figure('Position', normal_figure_size);
        % subplot(ceil(height(solver_tab)/2),2,ii)
        solver_name = solver_tab.solver_name(ii);
        hold on
        data_ii = [dtable.cpu_time_nlp(dtable.solver_name==solver_name,:)';dtable.cpu_time_lpec(dtable.solver_name==solver_name,:)'];
        % data_ii = log10(data_ii);
        bar(1:n_problems,data_ii,'stacked')
        legend_str = [append(solver_name,' NLP');append(solver_name,' LPEC')];
        legend(legend_str,'Location', 'northeast',BackgroundAlpha=0.5);
        xlabel('Problem instance')
        ylabel('Total wall time in seconds')
        ylim([0 cpu_max])
        % hold off
        figure(handle_fig_bar_timings);
        set(gca,'fontsize', FontSize);
        set(gca,'YScale','log');
        ax = gca;
        ax.ColorOrder= color_order;
        ax.LineStyleOrder = line_style_order;
        ax.LineStyleCyclingMethod = line_style_cycling_method;
        exportgraphics(gca, [filename '_bar_timing_' num2str(ii) '.pdf']);
    end
    % ;
end
%% Sucess faulire statistics
if settings.success_fail_statistics
    handle_fig_success_stats = figure('Position', normal_figure_size);
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        counts_ii = [sum(dtable.success(dtable.solver_name==solver_name,:));sum(dtable.max_iterations_reached(dtable.solver_name==solver_name,:));sum(dtable.problem_infeasible(dtable.solver_name==solver_name,:))];
        counts_ii = [counts_ii; sum(dtable.success(dtable.solver_name==solver_name,:)==0)-sum(counts_ii(2:3))];
        counts = [counts,counts_ii];
        legend_str = [legend_str ;solver_name];
    end
    figure(handle_fig_success_stats);
    ylabel('Count')
    hold on;
    bar(counts)
    xticks(1:4)
    xticklabels({'Sucess','Max. iters.','Infeasible','Other'});
    hold off
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northeast');
    if settings.save_plot
        exportgraphics(gca, [filename '_fail_stats.pdf']);
    end
end
%% Solved in phase i;
if settings.solved_in_phase_i
    % plot histogram of problems solved in presolve RELATIVE
    handle_fig_stat_point = figure('Position', normal_figure_size);
    ylabel('MPECs solved in Phase I [\%]')
    % ylabel('Count')
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        phase_i_success_count = dtable(dtable.solver_name==solver_name,:);
        counts = [counts,100*sum(phase_i_success_count.solved_in_phase_i)/length(phase_i_success_count.solved_in_phase_i)];
        legend_str = [legend_str ;solver_name];
        % solver_names = [solver_names;solver_name]'
    end
    figure(handle_fig_stat_point);
    hold on;
    bar(counts)
    xticks(1:height(solver_tab))
    xticklabels(legend_str)
    hold off
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northwest');
    exportgraphics(gca, [filename '_solved_in_phase_i_percent.pdf']);

    % plot histogram of problems solved in presolve ABSOLUTE
    handle_fig_stat_point_abs = figure('Position', normal_figure_size);
    ylabel('MPECs solved in Phase I')
    % ylabel('Count')
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        phase_i_success_count = dtable(dtable.solver_name==solver_name,:);
        counts = [counts,sum(phase_i_success_count.solved_in_phase_i)];
        legend_str = [legend_str ;solver_name];
        % solver_names = [solver_names;solver_name]'
    end
    figure(handle_fig_stat_point_abs);
    hold on;
    bar(counts)
    xticks(1:height(solver_tab))
    xticklabels(legend_str)
    hold off
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northwest');
    exportgraphics(gca, [filename '_solved_in_phase_i.pdf']);

end

%% Stationary points (todo only if sucessful)
if settings.stationary_points
    handle_fig_stat_point = figure('Position', normal_figure_size);
    xlabel('Stationary point')
    ylabel('Count')
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        stat_point_count = stationary_points_tab(stationary_points_tab.solver_name==solver_name,:);
        counts = [counts,stat_point_count.count];
        legend_str = [legend_str ;solver_name];
    end
    figure(handle_fig_stat_point);
    hold on;
    bar(counts)
    xticks(1:length(stationary_points_tab.multiplier_based_stationarity(stationary_points_tab.solver_name==solver_name,:)))
    xticklabels(stationary_points_tab.multiplier_based_stationarity(stationary_points_tab.solver_name==solver_name,:))
    hold off
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northwest');
    exportgraphics(gca, [filename '_stat_points.pdf']);
end
%% B-stationarity
if settings.b_stationarity
    handle_fig_b_stat = figure('Position', normal_figure_size);
    xlabel('B-stationarity')
    ylabel('Count')
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        stat_point_count = b_stationary_tab(b_stationary_tab.solver_name==solver_name,:);
        % stat_point_count = b_stationary_tab_among_sucess(b_stationary_tab_among_sucess.solver_name==solver_name,:);
        counts = [counts,stat_point_count.count];
        legend_str = [legend_str ;solver_name];
    end
    figure(handle_fig_b_stat);
    hold on;
    bar(counts)
    xticks(1:length(b_stationary_tab.b_stationarity(b_stationary_tab.solver_name==solver_name,:)))
    xticklabels(b_stationary_tab.b_stationarity(b_stationary_tab.solver_name==solver_name,:))
    hold off
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northwest');
    if settings.save_plot
        exportgraphics(gca, [filename '_B_stationarity.pdf']);
    end
end
%% B stationarity lpec residuals
if settings.b_stationarity
    handle_fig_b_lpec = figure('Position', wide_figure_size);
    ylabel('$|\nabla f(x^*)^\top d|$')
    xlabel('Problem instance')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_b_lpec);
        hold on;
        label = [char(solver_name)];
        obj_ii = abs(dtable.f_lpec(dtable.solver_name == solver_name));
        % if settings.plot_only_sucessful
        obj_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        obj_ii(obj_ii>0.99) = nan;
        obj_ii(obj_ii<1e-16) = 1e-16;
        % end
        stairs(obj_ii, 'LineWidth', linewidth, 'DisplayName', label);
    end
    problem_names = dtable.problem_name(dtable.solver_name == solver_name);
    hold off
    figure(handle_fig_b_lpec);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northwest',BackgroundAlpha=0.3);
    yline(1e-8,'k-','LineWidth', linewidth-0.5,'HandleVisibility','off')
    set(gca,'YScale','Log')
    ylim([0.5*1e-16 1])
    grid on
    if settings.save_plot
        exportgraphics(gca, [filename '_lpec_objective.pdf']);
    end
end

%% B stationary per solver
try

    if settings.b_stationarity
        plot_infeasible = 1;
        counts = [];
        legend_str = {};
        for ii=1:height(solver_tab)
            figure('Position', normal_figure_size);
            solver_name = solver_tab.solver_name(ii);
            stat_point_count = stationary_points_tab_with_b_stat(stationary_points_tab_with_b_stat.solver_name==solver_name,:);
            stat_point_count_is_not_b = stat_point_count.count(1:2:end);
            stat_point_count_is_b = stat_point_count.count(2:2:end);
            stat_point_count_is_not_b(end) = stat_point_count_is_not_b(end)+stat_point_count_is_b(end);
            stat_point_count_is_b(end) = 0;

            stat_point_count_is_not_b(end) = 0;

            x_ticks_length = length(stationary_points_tab_with_b_stat.multiplier_based_stationarity(stationary_points_tab_with_b_stat.solver_name==solver_name,:))/2;
            xtick_labels = stationary_points_tab.multiplier_based_stationarity(stationary_points_tab.solver_name==solver_name,:);
            hh = {};
            for jj = 1:length(xtick_labels)
                hh = [hh, xtick_labels(jj)];
            end
            hh{end}  = 'Failure';
            xtick_labels = hh;

            if ~plot_infeasible
                stat_point_count_is_not_b(end) = [];
                stat_point_count_is_b(end) = [];
                x_ticks_length = x_ticks_length-1;
            end
            bar([stat_point_count_is_b,stat_point_count_is_not_b],'stacked')
            title(solver_name)
            legend({'B stationary','not B stationary'},'location','northwest')
            xlabel('Stationary point')
            ylabel('Count')
            xticks(1:x_ticks_length)
            xticklabels(xtick_labels )
            set(gca,'fontsize', FontSize);
            % set(gca,'YScale','Log')
            if settings.save_plot
                exportgraphics(gca, [filename '_stat_points_details' num2str(ii) '.pdf']);
            end
        end
    end
catch
end

%% Objective function values

if settings.objective
    % plot objective values for all prob
    handle_fig_objective= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    % ylabel('Scaled objective value - $\mathrm{sign}(f) \log_{10}(|f|)$')
    ylabel('Objective value')
    % obj_cut_off = 10e3; % if larger then this in abs value, cap
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_objective);
        hold on;
        label = [char(solver_name)];
        obj_ii = dtable.f(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            obj_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        % obj_ii(obj_ii>obj_cut_off) = obj_cut_off;
        % obj_ii(obj_ii<-obj_cut_off) = -obj_cut_off;
        % obj_ii = sign(obj_ii ).*log10(abs(obj_ii )); % rescale to log
        stairs([obj_ii;obj_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    problem_names = dtable.problem_name(dtable.solver_name == solver_name);
    hold off
    figure(handle_fig_objective);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'southeast');
    if settings.save_plot
        exportgraphics(gca, [filename '_objective.pdf']);
    end
end


if settings.objective_rescaled
    % plot objective values for all prob
    handle_fig_objective= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Scaled objective value - $\mathrm{sign}(f) \log_{10}(|f|)$')
    % obj_cut_off = 10e3; % if larger then this in abs value, cap
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_fig_objective);
        hold on;
        label = [char(solver_name)];
        obj_ii = dtable.f(dtable.solver_name == solver_name);
        obj_ii = sign(obj_ii ).*log10(abs(obj_ii )); % rescale to log
        if settings.plot_only_sucessful
            obj_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([obj_ii;obj_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    problem_names = dtable.problem_name(dtable.solver_name == solver_name);
    hold off
    figure(handle_fig_objective);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'southeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_objective_rescaled.pdf']);
    end
end

%% Plot number of NLPs solves
if settings.nlps_solved
    handle_n_nlp_total= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('NLPs solved')
    if settings.bar_comparisson_plots
        nlps_solved = [];
        for ii=1:height(solver_tab)
            nlps_solved_ii = dtable.n_nlp_total(dtable.solver_name == solver_tab.solver_name(ii));
            if settings.plot_only_sucessful
                nlps_solved_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
            end
            nlps_solved  = [nlps_solved,nlps_solved_ii];
        end
        hold off;
        figure(handle_n_nlp_total);
        hold on;
        bar(nlps_solved);
        legend(solver_tab.solver_name)
    else
        for ii=1:height(solver_tab)
            solver_name = solver_tab.solver_name(ii);
            hold off;
            figure(handle_n_nlp_total);
            hold on;
            label = [char(solver_name)];
            n_nlp_total_ii = dtable.n_nlp_total(dtable.solver_name == solver_name);
            if settings.plot_only_sucessful
                n_nlp_total_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
            end
            stairs([n_nlp_total_ii;n_nlp_total_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
        hold off
    end
    figure(handle_n_nlp_total);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    max_nlp = max(dtable.n_nlp_total(dtable.success ==1))+2;
    ylim([0 max_nlp]);
    legend('Location', 'northeast','BackgroundAlpha',0.3,'NumColumns',2);
    if settings.save_plot
        exportgraphics(gca, [filename '_nlps_solved.pdf']);
    end
end

%% Plot max nlp wall time
if settings.nlp_cpu_time
    handle_nlp_cpu_time= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Aggregated NLP wall time')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_nlp_cpu_time);
        hold on;
        label = [char(solver_name)];
        nlp_cpu_time_ii = dtable.cpu_time_nlp(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            nlp_cpu_time_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([nlp_cpu_time_ii;nlp_cpu_time_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_nlp_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_nlp_cpu_time.pdf']);
    end
end


%% Plot  max nlp wall time
if settings.max_nlp_cpu_time
    handle_max_nlp_cpu_time= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Maximal NLPs wall time')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_max_nlp_cpu_time);
        hold on;
        label = [char(solver_name)];
        max_nlp_cpu_time_ii = dtable.max_cpu_time_nlp(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            max_nlp_cpu_time_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([max_nlp_cpu_time_ii;max_nlp_cpu_time_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_max_nlp_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')

    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_nlp_cpu_time.pdf']);
    end
end



%% Plot max nlp wall time phase i
if settings.max_nlp_cpu_time_phase_i
    handle_max_nlp_cpu_time= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Maximal NLPs wall time (Phase I)')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_max_nlp_cpu_time);
        hold on;
        label = [char(solver_name)];
        max_nlp_cpu_time_ii = dtable.max_cpu_time_nlp_phase_i(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            max_nlp_cpu_time_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([max_nlp_cpu_time_ii;max_nlp_cpu_time_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_max_nlp_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')

    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_nlp_cpu_time_phase_i.pdf']);
    end
end


%% Plot  max nlp wall time phase ii
if settings.max_nlp_cpu_time_phase_ii
    handle_max_nlp_cpu_time= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Maximal NLPs wall time (Phase II)')
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_max_nlp_cpu_time);
        hold on;
        label = [char(solver_name)];
        max_nlp_cpu_time_ii = dtable.max_cpu_time_nlp_phase_ii(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            max_nlp_cpu_time_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([max_nlp_cpu_time_ii;max_nlp_cpu_time_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_max_nlp_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')

    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_nlp_cpu_time_phase_ii.pdf']);
    end
end

%% Max LPEC wall time
if settings.max_lpec_cpu_time
    handle_max_lpec_cpu_time = figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Maximal LPECs wall time')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_max_lpec_cpu_time);
        hold on;
        label = [char(solver_name)];
        max_lpec_cpu_time_ii = dtable_lpec.max_cpu_time_lpec(dtable_lpec.solver_name == solver_name);
        if settings.plot_only_sucessful
            max_lpec_cpu_time_ii(dtable_lpec.success(dtable_lpec.solver_name == solver_name)==0) = nan;
        end
        stairs([max_lpec_cpu_time_ii;max_lpec_cpu_time_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_max_lpec_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')

    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'best',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_lpec_cpu_time.pdf']);
    end
end

% sorted
if settings.max_lpec_cpu_time
    handle_max_lpec_cpu_time = figure('Position', normal_figure_size);
    xlabel('Problem instance')
    ylabel('Maximal LPECs wall time')
    max_lpec_cpu_time = [];
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        max_lpec_cpu_time_ii = dtable_lpec.max_cpu_time_lpec(dtable_lpec.solver_name == solver_name);
        if settings.plot_only_sucessful
            max_lpec_cpu_time_ii(dtable_lpec.success(dtable_lpec.solver_name == solver_name)==0) = nan;
        end
        max_lpec_cpu_time = [max_lpec_cpu_time,max_lpec_cpu_time_ii];
    end
    % sort accodring to max
    [~,ii_max] = max(max(max_lpec_cpu_time));
    [~,ind_order] = sort(max_lpec_cpu_time(:,ii_max));
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_max_lpec_cpu_time);
        hold on;
        label = [char(solver_name)];
        stairs([max_lpec_cpu_time(ind_order,ii);max_lpec_cpu_time(ind_order(end),ii)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off

    figure(handle_max_lpec_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')

    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northwest',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_max_lpec_cpu_time_sorted.pdf']);
    end
end

%% Plot LPEC wall time
if settings.lpecs_cpu_time
    handle_lpecs_cpu_time= figure('Position', normal_figure_size);
    xlabel('Problem instance')
    ylabel('Total LPEC wall time')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_lpecs_cpu_time);
        hold on;
        label = [char(solver_name)];
        cpu_time_lpec_ii = dtable_lpec.cpu_time_lpec(dtable_lpec.solver_name == solver_name);
        cpu_time_lpec_ii  = cpu_time_lpec_ii +1e-6;
        if settings.plot_only_sucessful
            cpu_time_lpec_ii(dtable_lpec.success(dtable_lpec.solver_name == solver_name)==0) = nan;
        end
        stairs([cpu_time_lpec_ii;cpu_time_lpec_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end
    hold off
    figure(handle_lpecs_cpu_time);
    set(gca,'fontsize', FontSize);
    set(gca, 'YScale', 'log')
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northwest', BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_lpecs_cpu_time.pdf']);
    end
end

%% Plot LPEC vs NLP wall time

if settings.nlp_lpec_cpu_comparisson
    % for ii=1:height(solver_tab)
    %     % figure('Position', normal_figure_size);
    %     solver_name = solver_tab.solver_name(ii);
    %     hold on
    %     cpu_time_nlp_ii = success_tab.cpu_time_nlp(success_tab.solver_name==solver_name,:);
    %     cpu_time_lpec_ii = success_tab.cpu_time_lpec(success_tab.solver_name==solver_name,:);
    %     cpu_max = max([cpu_time_lpec_ii;cpu_time_nlp_ii]);
    %     cpu_min = max([min(cpu_time_lpec_ii),min(cpu_time_nlp_ii)]);
    %     tt = linspace(cpu_min*0.5,cpu_max*1.1,1e2);
    %     if 1
    %         figure('Position',normal_figure_size)
    %         loglog(tt,tt,'k','LineWidth',2)
    %         hold on
    %         loglog(cpu_time_nlp_ii,cpu_time_lpec_ii,'o','LineWidth',2)
    %         axis equal
    %         xlim([cpu_min*0.5 cpu_max*1.1])
    %         ylim([cpu_min*0.5 cpu_max*1.1])
    %
    %     else
    %         plot(cpu_time_nlp_ii,cpu_time_lpec_ii,'o','LineWidth',2)
    %         hold on
    %         plot(tt,tt,'k','LineWidth',2)
    %         % axis equal
    %         xlim([1e-6 cpu_max])
    %         ylim([1e-6 cpu_max])
    %     end
    %     ylabel('LPEC wall time')
    %     xlabel('NLP wall time')
    %     set(gca,'FontSize',FontSize);
    %     exportgraphics(gca, [filename '_nlp_vs_lpec_' num2str(ii) '.pdf']);
    % end

    % in one line
    f = figure('Position',wide_figure_size);
    for ii=1:height(solver_tab)
        subplot(1,height(solver_tab),ii)
        solver_name = solver_tab.solver_name(ii);
        title(solver_name )
        hold on
        cpu_time_nlp_ii = success_tab_lpec.cpu_time_nlp(success_tab_lpec.solver_name==solver_name,:);
        cpu_time_lpec_ii = success_tab_lpec.cpu_time_lpec(success_tab_lpec.solver_name==solver_name,:);

        ratio_cpu_lpec = sum(cpu_time_lpec_ii>cpu_time_nlp_ii)/length(cpu_time_lpec_ii)*100
        % fprintf(append(solver_name,": NLP faster than LPEC in %2.2f percent of cases \n"),ratio_cpu_lpec)
        fprintf(append(solver_name,": NLP faster than LPEC in %2.2f percent of cases \n"),ratio_cpu_lpec)
        cpu_max = max([cpu_time_lpec_ii;cpu_time_nlp_ii]);
        cpu_min = max([min(cpu_time_lpec_ii),min(cpu_time_nlp_ii)]);
        tt = linspace(cpu_min*0.5,cpu_max*1.1,1e2);
        if 1
            loglog(tt,tt,'k','LineWidth',2)
            hold on
            loglog(cpu_time_nlp_ii,cpu_time_lpec_ii,'.','LineWidth',2,'MarkerSize',10)
            axis equal
            xlim([cpu_min*0.5 cpu_max*1.1])
            ylim([cpu_min*0.5 cpu_max*1.1])
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
        else

            plot(tt,tt,'k','LineWidth',2)
            hold on
            plot(cpu_time_nlp_ii,cpu_time_lpec_ii,'.','LineWidth',2,'MarkerSize',10)
            % axis equal
            xlim([0 cpu_max])
            ylim([0 cpu_max])
        end
        if ii == 1
            ylabel('LPEC wall time')
        end
        xlabel('NLP wall time')
        set(gca,'FontSize',FontSize);
    end
    if settings.save_plot
        exportgraphics(f, [filename '_nlp_vs_lpec.pdf']);
    end
end

%% Plot number of LPECs solves
if settings.lpecs_solved
    handle_n_lpec_total = figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('LPECs solved')
    if settings.bar_comparisson_plots
        lpecs_solved = [];
        for ii=1:height(solver_tab)
            lpecs_solved_ii = dtable.n_lpec_total(dtable.solver_name == solver_tab.solver_name(ii));
            lpecs_solved  = [lpecs_solved,lpecs_solved_ii];
        end
        hold off;
        figure(handle_n_lpec_total);
        hold on;
        bar(lpecs_solved);
        legend(solver_tab.solver_name)
    else
        for ii=1:height(solver_tab)
            solver_name = solver_tab.solver_name(ii);
            hold off;
            figure(handle_n_lpec_total);
            hold on;
            label = [char(solver_name)];
            n_lpec_total_ii = dtable.n_lpec_total(dtable.solver_name == solver_name);
            if settings.plot_only_sucessful
                n_lpec_total_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
            end
            stairs([n_lpec_total_ii;n_lpec_total_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
    end
    problem_names = dtable.problem_name(dtable.solver_name == solver_name);
    hold off
    figure(handle_n_lpec_total);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    % xticks(1:length(problem_names));
    % xticklabels(problem_names);
    % xtickangle(45)
    % set(gca,xtickangle,45);
    if settings.save_plot
        exportgraphics(gca, [filename '_lpecs_solved.pdf']);
    end
end

%% Plot number of total active set changes during LPECs solves
if settings.active_set_changes
    handle_active_set_changes= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('Active set changes')

    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        hold off;
        figure(handle_active_set_changes);
        hold on;
        label = [char(solver_name)];
        active_set_changes_ii = dtable.active_set_changes(dtable.solver_name == solver_name);
        if settings.plot_only_sucessful
            active_set_changes_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
        end
        stairs([active_set_changes_ii;active_set_changes_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
    end

    problem_names = dtable.problem_name(dtable.solver_name == solver_name);
    hold off
    figure(handle_active_set_changes);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    % xticks(1:length(problem_names));
    % xticklabels(problem_names);
    % xtickangle(45)
    % set(gca,xtickangle,45);
    if settings.save_plot
        exportgraphics(gca, [filename '_active_set_changes.pdf']);
    end
end
%% Plot number of biactive constraints at solution;
if settings.n_biactive
    handle_n_biactive= figure('Position', wide_figure_size);
    xlabel('Problem instance')
    ylabel('$|\mathcal{I}_{00}|$')
    if settings.bar_comparisson_plots
        n_biactive = [];
        for ii=1:height(solver_tab)
            n_biactive_ii = dtable.n_biactive(dtable.solver_name == solver_tab.solver_name(ii));
            n_biactive  = [n_biactive,n_biactive_ii];
        end
        hold off;
        figure(handle_n_biactive);
        hold on;
        bar(n_biactive);
        legend(solver_tab.solver_name)
    else
        for ii=1:height(solver_tab)
            solver_name = solver_tab.solver_name(ii);
            hold off;
            figure(handle_n_biactive);
            hold on;
            label = [char(solver_name)];
            n_biactive_ii = dtable.n_biactive(dtable.solver_name == solver_name);
            if settings.plot_only_sucessful
                n_biactive_ii(dtable.success(dtable.solver_name == solver_name)==0) = nan;
            end
            stairs([n_biactive_ii;n_biactive_ii(end)], 'LineWidth', linewidth, 'DisplayName', label);
        end
    end
    hold off
    figure(handle_n_biactive);
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder= color_order;
    ax.LineStyleOrder = line_style_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    legend('Location', 'northeast',BackgroundAlpha=0.3);
    if settings.save_plot
        exportgraphics(gca, [filename '_biactive.pdf']);
    end
end

%% N biactive summary;
% if settings.n_biactive_count
if settings.n_biactive
    handle_fig_n_biactive_count = figure('Position', normal_figure_size);
    ylabel('MPECs with no biactive constraints [\%]')
    % ylabel('Count')
    counts = [];
    legend_str = {};
    for ii=1:height(solver_tab)
        solver_name = solver_tab.solver_name(ii);
        n_biactive_ii = dtable.n_biactive(dtable.solver_name == solver_tab.solver_name(ii));
        n_success_ii = sum(dtable.success(dtable.solver_name == solver_tab.solver_name(ii)));
        n_biactive_count = 100*sum(n_biactive_ii == 0)./n_success_ii;
        counts = [counts,n_biactive_count];
        legend_str = [legend_str ;solver_name];
        % solver_names = [solver_names;solver_name]'
    end
    figure(handle_fig_n_biactive_count);
    hold on;
    bar(counts)
    xticks(1:height(solver_tab))
    xticklabels(legend_str)
    hold off
    set(gca,'fontsize', FontSize);
    ax = gca;
    ax.ColorOrder = color_order;
    ax.LineStyleCyclingMethod = line_style_cycling_method;
    hold off;
    legend(legend_str,'Location', 'northwest');
    if settings.save_plot
        exportgraphics(gca, [filename '_n_biactive_count.pdf']);
    end
end

%%
if 0
    try
        solver_name_fast = solver_tab.solver_name(1);
        solver_name_slow = solver_tab.solver_name(2);
        cpu_time_fast = success_tab.cpu_time(success_tab.solver_name == solver_name_fast);
        cpu_time_slow  = success_tab.cpu_time(success_tab.solver_name == solver_name_slow);
        cpu_max = max([cpu_time_slow;cpu_time_fast]);
        tt = linspace(0,cpu_max,10);
        handle_compare= figure('Position', normal_figure_size);
        scatter(cpu_time_slow,cpu_time_fast,'LineWidth',2)
        hold on
        plot(tt,tt,'k','LineWidth',2)
        ylabel(solver_name_fast)
        xlabel(solver_name_slow)
        axis equal
        xlim([0 cpu_max])
        ylim([0 cpu_max])
        if settings.save_plot
            exportgraphics(gca, [filename '_compare.pdf']);
        end
    catch
    end
end
