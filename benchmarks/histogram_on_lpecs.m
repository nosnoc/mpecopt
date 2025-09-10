%% Evaluating LPEC MILP efficiency:
close all; clc;
latexify_plot()
% dtable;
% lpec_dstruct;

S = load('macmpec_general_01-Sep-2025_lpec_details');
lpec_dstruct = S.lpec_dstruct;
S = load('macmpec_general_01-Sep-2025');
dtable = S.dtable;


plot_cpu = false;
% Filter by solver and success

% solver_names  = ["Gurobi", "Gurobi-early-I", "Gurobi-early-I-II", "Highs", "Highs-early"];
solver_names = unique(dtable.solver_name);

solver = solver_names{1};
idx = (dtable.solver_name == solver) & dtable.success == 1;


% ------------------------
% Helper functions
% ------------------------

% Replace zeros by ones in a vector
function y = replace_zeros(x)
if all(x == 0)
    y = ones(size(x));
else
    x(x == 0) = 1;
    y = x;
end
end

% Postprocess nodecount cells
process_cells_cummulative = @(C) cellfun(@(x) ...
    (isempty(x) * 0) + ...             % empty → 0
    (~isempty(x) * (sum(replace_zeros(x)))), ...
    C);

process_cells = @(C) cellfun(@(x) ...
    (isempty(x) * 0) + ...             % empty → 0
    (~isempty(x) * ((replace_zeros(x)))), ...
    C, UniformOutput=false);

% Postprocess cpu time cells
process_cpu = @(C) cellfun(@(x) ...
    (isempty(x) * 0) + (~isempty(x) * sum(x)), ...
    C);

% Conversion of raw vector (per problem, for iteration plots)
function y = convert_vec(x)
if isempty(x)
    y = 0;
else
    if all(x == 0)
        y = ones(size(x));
    else
        x(x == 0) = 1;
        y = x;
    end
end
end

%%
%------------------------
% Nodecount statistics (for each lpec call seperatalyover solver calls)
% ------------------------
nc_i  = cell2mat(process_cells(lpec_dstruct.nodecount_phase_i(idx)));
nc_ii = cell2mat(process_cells(lpec_dstruct.nodecount_phase_ii(idx)));
nc_tot = [nc_i, nc_ii];

log_plot = true;
figure;

% --- Phase I ---
subplot(1,2,1)
if log_plot 
    h1 = histogram(log2(nc_i), 'BinMethod', 'integers');
    min_pow = floor(min(log2(nc_i)));
    max_pow = ceil(max(log2(nc_i)));
    set(gca, 'XTick', min_pow:max_pow)
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('$2^{%g}$', x), min_pow:max_pow, 'UniformOutput', false))
else
       h1 = histogram(nc_i, 'BinMethod', 'integers');
end

% Set ticks manually based on your data range
title('Nodecount Phase I - each LPEC call');
xlabel('Nodes'); ylabel('Frequency'); grid on;
xticks = get(gca, 'XTick');
% annotate counts
for k = 1:numel(h1.BinEdges)-1
    if h1.Values(k) > 0
        x = (h1.BinEdges(k) + h1.BinEdges(k+1))/2;
        y = h1.Values(k);
        text(x, y, num2str(y), 'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end
% adjust ylim
ymax1 = max(h1.Values);
ylim([0 ymax1 + max(2,ceil(0.2*ymax1))])

% --- Phase II ---
subplot(1,2,2)
if log_plot 
    h2 = histogram(log2(nc_ii), 'BinMethod', 'integers');
    % min_pow = floor(min(log2(nc_ii)));
    % max_pow = ceil(max(log2(nc_ii)));
    set(gca, 'XTick', min_pow:max_pow)
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('$2^{%g}$', x), min_pow:max_pow, 'UniformOutput', false))
    % enforce nicer axis if only one bar
    if numel(unique(nc_ii)) == 1
            xlim([-1 max_pow]);
    end
else
    h2 = histogram(nc_ii, 'BinMethod', 'integers');
    % enforce nicer axis if only one bar
if numel(unique(nc_ii)) == 1
    xlim([0 10])
end
end
title('Nodecount Phase II - each LPEC call');
xlabel('Nodes'); ylabel('Frequency'); grid on;


% annotate counts
for k = 1:numel(h2.BinEdges)-1
    if h2.Values(k) > 0
        x = (h2.BinEdges(k) + h2.BinEdges(k+1))/2;
        y = h2.Values(k);
        text(x, y, num2str(y), 'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

% adjust ylim
ymax2 = max(h2.Values);
ylim([0 ymax2 + max(2,ceil(0.2*ymax2))])

%%
% xticks(min(nc_ii):max(nc_ii));

% subplot(133)
% histogram(nc_tot, 'BinMethod', 'integers');
% title('Total Nodecount - each LPEC call'); xlabel('Nodes'); ylabel('Frequency'); grid on;
% xticks(min(nc_tot):max(nc_tot));


% ------------------------
% Nodecount statistics (cummulative over solver calls)
% ------------------------
nc_i  = process_cells_cummulative(lpec_dstruct.nodecount_phase_i(idx));
nc_ii = process_cells_cummulative(lpec_dstruct.nodecount_phase_ii(idx));
nc_tot = nc_i + nc_ii;
if 0
    % Histograms (nodecounts)
    figure;
    subplot(131)
    histogram(nc_i, 'BinMethod', 'integers');
    title('Nodecount Phase I'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_i):max(nc_i));

    subplot(132)
    histogram(nc_ii, 'BinMethod', 'integers');
    title('Nodecount Phase II'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_ii):max(nc_ii));

    subplot(133)
    histogram(nc_tot, 'BinMethod', 'integers');
    title('Total Nodecount'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_tot):max(nc_tot));
end
% Largest problem by total nodecount
[sort_val, sort_idx_rel] = sort(nc_tot);
max_idx_rel = sort_idx_rel(end-1);
% [~, max_idx_rel] = max(nc_tot);
max_idx_all = find(idx); % indices in original table
max_idx = max_idx_all(max_idx_rel);

fprintf('Problem with largest total nodecount is %s with %d nodes.\n', ...
    dtable.problem_name{max_idx}, nc_tot(max_idx_rel));

% ------------------------
% CPU time statistics
% ------------------------
cpu_i  = process_cpu(lpec_dstruct.cpu_time_lpec_phase_i(idx));
cpu_ii = process_cpu(lpec_dstruct.cpu_time_lpec_phase_ii(idx));
cpu_tot = cpu_i + cpu_ii;

% Histograms (CPU times)
if plot_cpu
    figure;
    subplot(131)
    histogram(cpu_i);
    title('CPU Time Phase I'); xlabel('Time [s]'); ylabel('Frequency'); grid on;

    subplot(132)
    histogram(cpu_ii);
    title('CPU Time Phase II'); xlabel('Time [s]'); ylabel('Frequency'); grid on;

    subplot(133)
    histogram(cpu_tot);
    title('Total CPU Time'); xlabel('Time [s]'); ylabel('Frequency'); grid on;
end

% ------------------------
% max_idx = 14;
%% Iteration plot for max problem (nodecounts)
% ------------------------
vec_i =  cell2mat(process_cells(lpec_dstruct.nodecount_phase_i(max_idx)));
vec_ii = cell2mat(process_cells(lpec_dstruct.nodecount_phase_ii(max_idx)));
% vec_i = convert_vec(vec_i);
% vec_ii = convert_vec(vec_ii);
x_i = 1:length(vec_i);
x_ii = (length(vec_i)+1):(length(vec_i)+length(vec_ii));
figure;
subplot(121)
plot(x_i, vec_i, '-o', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName','Phase I'); hold on;
plot(x_ii, vec_ii, '-s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName','Phase II');
% Connect the two phases
if ~isempty(vec_i) && ~isempty(vec_ii)
    plot([x_i(end), x_ii(1)], [vec_i(end), vec_ii(1)], '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
end
ylim([0 (max([vec_i, vec_ii])+1)*1.1])
xlabel('Iteration');
ylabel('Nodecount');
% title(sprintf('Nodecount per Iteration (%s)', dtable.problem_name{max_idx}));
legend('show'); grid on;
% ------------------------
% Iteration plot for max problem (CPU times)
% ------------------------
cpu_vec_i = lpec_dstruct.cpu_time_lpec_phase_i{max_idx};
cpu_vec_ii = lpec_dstruct.cpu_time_lpec_phase_ii{max_idx};
if isempty(cpu_vec_i), cpu_vec_i = 0; end
if isempty(cpu_vec_ii), cpu_vec_ii = 0; end
x_ci = 1:length(cpu_vec_i);
x_cii = (length(cpu_vec_i)+1):(length(cpu_vec_i)+length(cpu_vec_ii));
subplot(122)
plot(x_ci, cpu_vec_i, '-o', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName','Phase I'); hold on;
plot(x_cii, cpu_vec_ii, '-s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName','Phase II');
% Connect the two phases
if ~isempty(cpu_vec_i) && ~isempty(cpu_vec_ii)
    plot([x_ci(end), x_cii(1)], [cpu_vec_i(end), cpu_vec_ii(1)], '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
end
xlabel('Iteration');
ylabel('CPU Time [s]');
legend('show'); grid on;
%%
%% Iteration plot for max problem (nodecounts)
% ------------------------
vec_i = lpec_dstruct.nodecount_phase_i{max_idx};
vec_ii = lpec_dstruct.nodecount_phase_ii{max_idx};
vec_i = convert_vec(vec_i);
vec_ii = convert_vec(vec_ii);

% Combine vectors and create x-axis
combined_vec = [vec_i, vec_ii];
x_combined = 1:length(combined_vec);

% Create color array - blue for Phase I, orange for Phase II
colors = [repmat([0 0.4470 0.7410], length(vec_i), 1); ...
    repmat([0.8500 0.3250 0.0980], length(vec_ii), 1)];
%%
figure;
subplot(121)
b = bar(x_combined, combined_vec, 'FaceColor', 'flat');
b.CData = colors;
ylim([0 (max(combined_vec)+1)*1.1])
xlabel('Iteration');
ylabel('Nodecount');
% title(sprintf('Nodecount per Iteration (%s)', dtable.problem_name{max_idx}));
grid on;

% Add legend manually
hold on;
h1 = bar(NaN, NaN, 'FaceColor', [0 0.4470 0.7410]);
h2 = bar(NaN, NaN, 'FaceColor', [0.8500 0.3250 0.0980]);
legend([h1, h2], {'Phase I', 'Phase II'}, 'Location', 'best');

% ------------------------
% Iteration plot for max problem (CPU times)
% ------------------------
cpu_vec_i = lpec_dstruct.cpu_time_lpec_phase_i{max_idx};
cpu_vec_ii = lpec_dstruct.cpu_time_lpec_phase_ii{max_idx};

cpu_vec_i(cpu_vec_i==0)  = [];
cpu_vec_ii(cpu_vec_ii==0)  = [];

if isempty(cpu_vec_i), cpu_vec_i = 0; end
if isempty(cpu_vec_ii), cpu_vec_ii = 0; end

% Ensure CPU vectors are column vectors and combine them
cpu_vec_i = cpu_vec_i(:);  % Convert to column vector
cpu_vec_ii = cpu_vec_ii(:);  % Convert to column vector
combined_cpu = [cpu_vec_i; cpu_vec_ii];
x_cpu_combined = 1:length(combined_cpu);

% Create color array for CPU plot
cpu_colors = [repmat([0 0.4470 0.7410], length(cpu_vec_i), 1); ...
    repmat([0.8500 0.3250 0.0980], length(cpu_vec_ii), 1)];

subplot(122)
b2 = bar(x_cpu_combined, combined_cpu, 'FaceColor', 'flat');
b2.CData = cpu_colors;
xlabel('Iteration');
ylabel('CPU Time [s]');
grid on;

% Add legend manually
hold on;
h3 = bar(NaN, NaN, 'FaceColor', [0 0.4470 0.7410]);
h4 = bar(NaN, NaN, 'FaceColor', [0.8500 0.3250 0.0980]);
legend([h3, h4], {'Phase I', 'Phase II'}, 'Location', 'best');

