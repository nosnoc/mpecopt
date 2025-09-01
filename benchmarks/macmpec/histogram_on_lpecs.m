%% Evaluating LPEC MILP efficiency:
close all; clc;
dtable;
lpec_dstruct;

plot_cpu = false;
% Filter by solver and success
solver = "MPECopt-Reg-Early";
solver = "MPECopt-Reg-Gurobi";
idx = strcmp(dtable.solver_name, solver) & dtable.success == 1;

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

% ------------------------
% Nodecount statistics (for each lpec call seperatalyover solver calls)
% ------------------------
nc_i  = cell2mat(process_cells(lpec_dstruct.nodecount_phase_i(idx)));
nc_ii = cell2mat(process_cells(lpec_dstruct.nodecount_phase_ii(idx)));
nc_tot = [nc_i, nc_ii];

% Histograms (nodecounts)
figure; 
subplot(131)
histogram(nc_i, 'BinMethod', 'integers');
title('Nodecount Phase I - each LPEC call'); xlabel('Nodes'); ylabel('Frequency'); grid on;
% xticks(min(nc_i):max(nc_i));

subplot(132)
histogram(nc_ii, 'BinMethod', 'integers');
title('Nodecount Phase II - each LPEC call'); xlabel('Nodes'); ylabel('Frequency'); grid on;
xticks(min(nc_ii):max(nc_ii));

subplot(133)
histogram(nc_tot, 'BinMethod', 'integers');
title('Total Nodecount - each LPEC call'); xlabel('Nodes'); ylabel('Frequency'); grid on;
% xticks(min(nc_tot):max(nc_tot));



% ------------------------
% Nodecount statistics (cummulative over solver calls)
% ------------------------
nc_i  = process_cells_cummulative(lpec_dstruct.nodecount_phase_i(idx));
nc_ii = process_cells_cummulative(lpec_dstruct.nodecount_phase_ii(idx));
nc_tot = nc_i + nc_ii;

% Histograms (nodecounts)
figure; 
subplot(131)
histogram(nc_i, 'BinMethod', 'integers');
title('Nodecount Phase I'); xlabel('Nodes'); ylabel('Frequency'); grid on;
% xticks(min(nc_i):max(nc_i));

subplot(132)
histogram(nc_ii, 'BinMethod', 'integers');
title('Nodecount Phase II'); xlabel('Nodes'); ylabel('Frequency'); grid on;
xticks(min(nc_ii):max(nc_ii));

subplot(133)
histogram(nc_tot, 'BinMethod', 'integers');
title('Total Nodecount'); xlabel('Nodes'); ylabel('Frequency'); grid on;
% xticks(min(nc_tot):max(nc_tot));

% Largest problem by total nodecount
[~, max_idx_rel] = max(nc_tot);
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
%% Iteration plot for max problem (nodecounts)
% ------------------------
vec_i  = lpec_dstruct.nodecount_phase_i{max_idx};
vec_ii = lpec_dstruct.nodecount_phase_ii{max_idx};

vec_i  = convert_vec(vec_i);
vec_ii = convert_vec(vec_ii);

x_i  = 1:length(vec_i);
x_ii = (length(vec_i)+1):(length(vec_i)+length(vec_ii));

figure;
subplot(121)
stairs(x_i, vec_i, '-o', 'LineWidth', 1.5, 'DisplayName','Phase I'); hold on;
stairs(x_ii, vec_ii, '-s', 'LineWidth', 1.5, 'DisplayName','Phase II');
ylim([0 (max([vec_i, vec_ii])+1)*1.1])
xlabel('Iteration');
ylabel('Nodecount');
% title(sprintf('Nodecount per Iteration (%s)', dtable.problem_name{max_idx}));
legend('show'); grid on;

% ------------------------
% Iteration plot for max problem (CPU times)
% ------------------------
cpu_vec_i  = lpec_dstruct.cpu_time_lpec_phase_i{max_idx};
cpu_vec_ii = lpec_dstruct.cpu_time_lpec_phase_ii{max_idx};

if isempty(cpu_vec_i), cpu_vec_i = 0; end
if isempty(cpu_vec_ii), cpu_vec_ii = 0; end

x_ci  = 1:length(cpu_vec_i);
x_cii = (length(cpu_vec_i)+1):(length(cpu_vec_i)+length(cpu_vec_ii));

subplot(122)
stairs(x_ci, cpu_vec_i, '-o', 'LineWidth', 1.5, 'DisplayName','Phase I'); hold on;
stairs(x_cii, cpu_vec_ii, '-s', 'LineWidth', 1.5, 'DisplayName','Phase II');
xlabel('Iteration');
ylabel('CPU Time [s]');
title(sprintf('CPU Time per Iteration (%s)', dtable.problem_name{max_idx}));
legend('show'); grid on;
