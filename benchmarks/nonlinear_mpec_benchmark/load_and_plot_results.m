close all;

plot_settings.relative = 1;
plot_settings.relative_phase_i = 0;
plot_settings.relative_phase_ii = 0;
plot_settings.relative_lpec = 0;

plot_settings.absolute = 1;
plot_settings.absolute_phase_i = 0;
plot_settings.absolute_phase_ii = 0;
plot_settings.absolute_lpec = 0;

plot_settings.success_fail_statistics = 1;

plot_settings.nlp_lpec_cpu_comparisson = 0;

plot_settings.lpecs_solved = 1;
plot_settings.nlps_solved = 1;

plot_settings.lpecs_cpu_time = 0;
plot_settings.nlp_cpu_time = 0; % aggegated nlp times phase I and II

plot_settings.max_nlp_cpu_time = 1;
plot_settings.max_lpec_cpu_time = 0;

plot_settings.n_biactive_count = 0;
plot_settings.n_biactive = 1;
plot_settings.active_set_changes = 0;
plot_settings.bar_timing_plots = 0;
plot_settings.bar_comparisson_plots = 0;

plot_settings.max_nlp_cpu_time_phase_i = 0;
plot_settings.max_nlp_cpu_time_phase_ii = 0;
plot_settings.nlp_phases_cpu_time = 0;
plot_settings.lpec_phases_cpu_time = 0;
plot_settings.solved_in_phase_i = 0;

plot_settings.objective = 0;
plot_settings.objective_rescaled = 1;

% split cpu time among phases and show in bar plot (just report for mpecopt
% guroi)
plot_settings.stationary_points = 0;
plot_settings.b_stationarity = 0;
plot_settings.b_stationarty_as_success_criterion = 0;
plot_settings.plot_only_sucessful = 0;

plot_settings.save_plot = 1;

%%
results = 1;
switch results
    case 1
        S = load('nonlinear_mpec_large1_15-Sep-2025_5');
        % S = load('nonlinear_mpec_general2_09-Nov-2024');
        try
            dtable = S.dtable;
        catch
            dtable = S.dtable1;
        end
        solver_names  = ["MPECopt-Reg-Gurobi", "MPECopt-Reg-Gurobi-ET","Reg", "NLP", "$\ell_1$-Penalty"];

        filename = 'nonlinear_mpec';
        N_plot = [1:5];
        plot_settings.stationary_points = 0;
        plot_settings.b_stationarity = 0;
        plot_settings.n_biactive = 1;
        %% Who is not B stationary?
        if 0
            dtable_ii = dtable(dtable.solver_name == solver_names{1},:);
            dtable_jj = dtable(dtable.solver_name == solver_names{4},:);
            ind_not_B = find(dtable_jj.success == 1 & dtable_jj.b_stationarity == 0);
            % For Reg
            % diff in obj
            delta_f = dtable_ii.f(ind_not_B)-dtable_jj.f(ind_not_B);
            delta_f_relative_small = 100*abs(delta_f)./abs(dtable_ii.f(ind_not_B)+1e-16)
            delta_f_relative_small_ind = abs(delta_f)./abs(dtable_ii.f(ind_not_B)+1e-16)<1e-3;
            dtable_temp = dtable_jj;
            dtable_temp.b_stationarity(ind_not_B(delta_f_relative_small_ind)) = 1;
            dtable_temp.multiplier_based_stationarity(ind_not_B(delta_f_relative_small_ind)) = dtable_ii.multiplier_based_stationarity(ind_not_B(delta_f_relative_small_ind));
            dtable_temp.problem_name(ind_not_B(~delta_f_relative_small_ind))
            dtable_temp.prob_num(ind_not_B(~delta_f_relative_small_ind))
            dtable_temp.f_lpec(ind_not_B(~delta_f_relative_small_ind))
            dtable(dtable.solver_name == solver_names{4},:) = dtable_temp;
            % For pen
            dtable_ii = dtable(dtable.solver_name == solver_names{3},:);
            dtable_jj = dtable(dtable.solver_name == solver_names{6},:);
            ind_not_B = find(dtable_jj.success == 1 & dtable_jj.b_stationarity == 0);
            % diff in obj
            delta_f = dtable_ii.f(ind_not_B)-dtable_jj.f(ind_not_B);
            delta_f_relative_small = 100*abs(delta_f)./abs(dtable_ii.f(ind_not_B)+1e-16)
            delta_f_relative_small_ind = abs(dtable_ii.f(ind_not_B)-dtable_jj.f(ind_not_B))./abs(dtable_ii.f(ind_not_B)+1e-16)<1e-3;
            dtable_temp = dtable_jj;
            dtable_temp.b_stationarity(ind_not_B(delta_f_relative_small_ind)) = 1;
            dtable_temp.multiplier_based_stationarity(ind_not_B(delta_f_relative_small_ind)) = dtable_ii.multiplier_based_stationarity(ind_not_B(delta_f_relative_small_ind));
            dtable_temp.problem_name(ind_not_B(~delta_f_relative_small_ind))
            dtable_temp.prob_num(ind_not_B(~delta_f_relative_small_ind))
            dtable_temp.f_lpec(ind_not_B(~delta_f_relative_small_ind))
            dtable(dtable.solver_name == solver_names{6},:) = dtable_temp;
        end
    case 2

        S = load('macmpec_phase_i_29-Oct-2024');
        dtable = S.dtable;
        solver_names = ["Reg-LPEC", "Reg-Simple",...
            "Pen-$\ell_{\infty}$-LPEC","Pen-$\ell_{1}$-LPEC",...
            "Feasibility-$\ell_1$", "Feasibility-$\ell_{\infty}$"];
        filename = 'macmpec_phase_i';
        N_plot = [1,2,3,5,6];
        plot_settings.relative_phase_i = 1;
        plot_settings.absolute_phase_i = 1;

    case 3
        % todo: compine with highs and gurobi from general;
        S = load('macmpec_lpec_08-Nov-2024_3');
        try
            dtable = S.dtable;
        catch
            dtable = S.dtable1;
        end
        solver_names = ["Gurobi-MILP", "HiGHS-MILP",...
            "Reg-MPEC","$\ell_{1}$-MPEC",...
            "$\ell_{\infty}$-MPEC"];
        N_plot = 1:4;
        % N_plot = [3,5];
        plot_settings.relative= 0;
        plot_settings.absolute = 0;
        plot_settings.relative_lpec = 1;
        plot_settings.absolute_lpec = 1;
        plot_settings.nlp_lpec_cpu_comparisson = 1;
        plot_settings.max_lpec_cpu_time = 1;
    case 4
        S = load('debug5_30-Oct-2024.mat');
        dtable = S.dtable;
        filename = 'macmpec';

        solver_names  = ["MPECopt-Gurobi-def", "gurobi-not-tigher", "gurobi-phI-1e-2",...
            "gurobi not tight 1e-2", "highs-1e-3","highs-1e-3-dont try","reg-no recovery","reg-default",...
            'Reg1e-9','Reg1e-10'];

        N_plot = [1,2];

end

% Modift possibly which solvers are plotted
[~,temp_b] = unique(dtable.solver_name);
solver_name = dtable.solver_name(sort(temp_b));
% solver_tab = table(solver_name);

dtable_plot = [];
for ii = N_plot
    dtable_plot  = [dtable_plot; dtable(dtable.solver_name == solver_names(ii),:)];
end

%% Export to excel
write_in_excel_table = false;
if write_in_excel_table
    write_in_separate_sheets = true;
    filename_excel = [filename '_phase_i_details.xlsx'];
    if write_in_separate_sheets
        solver_tab = pivot(dtable_plot, Rows= "solver_name");
        for ii = 1:length(solver_tab.solver_name)
            solver_name =solver_tab.solver_name{ii};
            dtable_ii = dtable_plot(dtable_plot.solver_name == solver_name,:);
            writetable(dtable_ii,filename_excel,'Sheet',ii)
        end
    else
        writetable(dtable_plot,filename_excel,'Sheet',1)
    end
end
%% Plot results
plot_mpec_benchmark_result(dtable_plot,filename,plot_settings)

%%
dtable_jj = dtable_plot(dtable_plot.solver_name == solver_names{2},:);
hh = dtable_jj.success == 0;
dtable_jj.problem_name(find(hh));
% dtable_jj.problem_name(find(hh))


dtable_ii = dtable_plot(dtable_plot.solver_name == solver_names{5},:);
hh = dtable_ii.b_stationarity == 0 & dtable_ii.multiplier_based_stationarity == 'S';
dtable_ii.problem_name(find(hh));
find(hh')
