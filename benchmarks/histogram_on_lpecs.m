%% Evaluating LPEC MILP efficiency:
close all; clc;
latexify_plot()
% dtable;
% lpec_dstruct;

S = load('random_lpecs2_07-Sep-2025_lpec_details');
lpec_dstruct = S.lpec_dstruct;
S = load('random_lpecs2_07-Sep-2025');
%
% S = load('nonlinear_mpec_med_15-Sep-2025_lpec_details');
% lpec_dstruct = S.lpec_dstruct;
% S = load('nonlinear_mpec_med_15-Sep-2025');


dtable = S.dtable;

plot_cpu = false;
% Filter by solver and success

% solver_names  = ["Gurobi", "Gurobi-early-I", "Gurobi-early-I-II", "Highs", "Highs-early"];
solver_names = unique(dtable.solver_name);

solver = solver_names{3}
idx = (dtable.solver_name == solver) & dtable.success == 1;

fields = fieldnames(lpec_dstruct);

for i = 1:length(fields)
    lpec_dstruct_filtered.(fields{i}) = lpec_dstruct.(fields{i})(idx);
end

lpec_dstruct = lpec_dstruct_filtered;



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
nc_i  = cell2mat(process_cells(lpec_dstruct.nodecount_phase_i));
nc_ii = cell2mat(process_cells(lpec_dstruct.nodecount_phase_ii));
nc_tot = [nc_i, nc_ii];

log_plot = true;
figure;

% --- Phase I ---
subplot(1,2,1)
if log_plot
    h1 = histogram(log2(nc_i), 'BinMethod', 'integers');
    min_pow = floor(min(log2(nc_i)));
    max_pow = max(5,ceil(max(log2(nc_i))));
    set(gca, 'XTick', min_pow:max_pow)
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('$2^{%g}$', x), min_pow:max_pow, 'UniformOutput', false))
    if numel(unique(nc_ii)) == 1
        xlim([-1 max_pow]);
    end
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
% xlim([0 100])

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
nc_cummulative_i  = process_cells_cummulative(lpec_dstruct.nodecount_phase_i);
nc_cummulative_ii = process_cells_cummulative(lpec_dstruct.nodecount_phase_ii);
nc_tot_cummulative = nc_cummulative_i + nc_cummulative_ii;
if 0
    % Histograms (nodecounts)
    figure;
    subplot(131)
    histogram(nc_cummulative_i, 'BinMethod', 'integers');
    title('Nodecount Phase I'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_i):max(nc_i));

    subplot(132)
    histogram(nc_cummulative_ii, 'BinMethod', 'integers');
    title('Nodecount Phase II'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_ii):max(nc_ii));

    subplot(133)
    histogram(nc_tot_cummulative, 'BinMethod', 'integers');
    title('Total Nodecount'); xlabel('Nodes'); ylabel('Frequency'); grid on;
    % xticks(min(nc_tot):max(nc_tot));
end
%% Largest problem by total nodecount
[sort_val, sort_idx_rel] = sort(nc_tot_cummulative);
max_idx_rel = sort_idx_rel(end-2);
% [~, max_idx_rel] = max(nc_tot);
% max_idx_all = find(idx); % indices in original table
% max_idx = max_idx_all(max_idx_rel);
max_idx = max_idx_rel;
% max_idx = 167;

% [max_idx max_idx_rel]

% pack-rig1-8.nl
% max_idx = 128;
% max_idx_rel = 128;

% qpec-200-3
% max_idx = 166;
% max_idx_rel = 157;

% qpec-200-4
% max_idx = 167;
% max_idx_rel = 158;



%
% max_idx = 127;
% max_idx = 158;
% max_idx = 167;

fprintf('Problem with largest total nodecount is %s with %d nodes.\n', ...
    dtable.problem_name{max_idx}, nc_tot_cummulative(max_idx));

% ------------------------
% CPU time statistics
% ------------------------
cpu_i  = process_cpu(lpec_dstruct.cpu_time_lpec_phase_i);
cpu_ii = process_cpu(lpec_dstruct.cpu_time_lpec_phase_ii);
cpu_tot = cpu_i + cpu_ii;


% ------------------------
%% Iteration plot for max problem (nodecounts)
% ------------------------
vec_i =  cell2mat(process_cells(lpec_dstruct.nodecount_phase_i(max_idx)));
vec_ii = cell2mat(process_cells(lpec_dstruct.nodecount_phase_ii(max_idx)));
% vec_i = convert_vec(vec_i);
% vec_ii = convert_vec(vec_ii);
x_i = 1:length(vec_i);
x_ii = (length(vec_i)+1):(length(vec_i)+length(vec_ii));

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


%%

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