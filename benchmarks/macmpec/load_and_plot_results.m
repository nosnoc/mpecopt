close all;

plot_settings.relative = 1;
plot_settings.relative_phase_i = 0;
plot_settings.relative_phase_ii = 0;
plot_settings.relative_lpec = 0;

plot_settings.absolute = 1;
plot_settings.absolute_phase_i = 0;
plot_settings.absolute_phase_ii = 0;
plot_settings.absolute_lpec = 0;

plot_settings.success_fail_statistics = 0;
plot_settings.nlp_lpec_cpu_comparisson = 0;

plot_settings.lpecs_solved = 0;
plot_settings.nlps_solved = 0;

plot_settings.lpecs_cpu_time = 0;
plot_settings.nlp_cpu_time = 0; % aggegated nlp times phase I and II

plot_settings.max_nlp_cpu_time = 0;
plot_settings.max_lpec_cpu_time = 0;

plot_settings.n_biactive_count = 1;
plot_settings.n_biactive = 0;
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
plot_settings.plot_only_sucessful = 1;

plot_settings.save_plot = 1;

%%
results = 1;
switch results
    case 1
        % S = load('macmpec_general_30-Oct-2024');      
        % S = load('macmpec_general_07-Nov-2024');
        % S = load('macmpec_general_14-Sep-2025');
        S = load('macmpec_general_22-Sep-2025');
        dtable = S.dtable;
        solver_names = unique(dtable.solver_name);
        N_plot = [4 3 1 7 6 2];

        % solver_names  = ["MPECopt-Reg-Gurobi", "MPECopt-$\ell_1$-Gurobi", "MPECopt-Reg-Guroby-ET", ...
        %           "Reg", "NLP", ...
        %           "MINLP"];


        filename = 'macmpec';
        % N_plot = [1:6];
        % plot_settings.stationary_points = 1;
        % plot_settings.b_stationarity = 1;
        % plot_settings.n_biactive = 1;
        % plot_settings.solved_in_phase_i = 1;
     %% Post process to verify success and B-stationarity
     if 1
         ind_ref = 1;
         ind_fix = 4;
         dtable_ii = dtable(dtable.solver_name == solver_names{ind_ref},:); % Reference solution
         dtable_jj = dtable(dtable.solver_name == solver_names{ind_fix},:); % To be chechked
         ind_not_success = find(dtable_ii.success == 1 & dtable_jj.success == 0);
         ind_not_success = find(dtable_jj.success == 0);
         % For Reg
         % diff in obj
         delta_f = dtable_ii.f(ind_not_success)-dtable_jj.f(ind_not_success);
         delta_f_relative_small = 100*abs(delta_f)./abs(dtable_ii.f(ind_not_success)+1e-16);
         delta_f_relative_small_ind = abs(delta_f)./abs(dtable_ii.f(ind_not_success)+1e-16)<1e-2;
         dtable_temp = dtable_jj;
         dtable_temp.success(ind_not_success(delta_f_relative_small_ind)) = dtable_ii.success(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.b_stationarity(ind_not_success(delta_f_relative_small_ind)) =   dtable_ii.b_stationarity(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.multiplier_based_stationarity(ind_not_success(delta_f_relative_small_ind)) = dtable_ii.multiplier_based_stationarity(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.problem_name(ind_not_success(~delta_f_relative_small_ind))
         dtable_temp.prob_num(ind_not_success(~delta_f_relative_small_ind))
         dtable_temp.f_lpec(ind_not_success(~delta_f_relative_small_ind))
         dtable_jj = dtable_temp;
         dtable(dtable.solver_name == solver_names{ind_fix},:) = dtable_temp;

         % still not successful
         ind_not_success = find(dtable_ii.success == 1 & dtable_jj.success == 0);
         dtable_jj(ind_not_success,:)

         % For pen
         ind_ref = 3;
         ind_fix = 1;
         dtable_ii = dtable(dtable.solver_name == solver_names{ind_ref},:); % Reference solution
         dtable_jj = dtable(dtable.solver_name == solver_names{ind_fix},:); % To be chechked
         ind_not_success = find(dtable_ii.success == 1 & dtable_jj.success == 0);
         % diff in obj
         delta_f = dtable_ii.f(ind_not_success)-dtable_jj.f(ind_not_success);
         delta_f_relative_small = 100*abs(delta_f)./abs(dtable_ii.f(ind_not_success)+1e-16);
         delta_f_relative_small_ind = abs(dtable_ii.f(ind_not_success)-dtable_jj.f(ind_not_success))./abs(dtable_ii.f(ind_not_success)+1e-16)<1e-2;
         dtable_temp = dtable_jj;
         dtable_temp.success(ind_not_success(delta_f_relative_small_ind)) = dtable_ii.success(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.b_stationarity(ind_not_success(delta_f_relative_small_ind)) =   dtable_ii.b_stationarity(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.multiplier_based_stationarity(ind_not_success(delta_f_relative_small_ind)) = dtable_ii.multiplier_based_stationarity(ind_not_success(delta_f_relative_small_ind));
         dtable_temp.problem_name(ind_not_success(~delta_f_relative_small_ind))
         dtable_temp.prob_num(ind_not_success(~delta_f_relative_small_ind))
         dtable_temp.f_lpec(ind_not_success(~delta_f_relative_small_ind))
         dtable_jj = dtable_temp;
         dtable(dtable.solver_name == solver_names{ind_fix},:) = dtable_temp;

         % still not successful
         ind_not_success = find(dtable_ii.success == 1 & dtable_jj.success == 0);
         dtable_jj(ind_not_success,:)

         % number of feasible but not b stat in minlp
         ind_minlp = 2;
         dtable_ii = dtable(dtable.solver_name == solver_names{ind_minlp},:);
         sum(dtable_ii.success == 0 & dtable_ii.problem_infeasible == 0)
       
     end


    case 2

        S = load('macmpec_phase_i_08-Nov-2024');
        dtable = S.dtable;
        solver_names = ["Reg-LPEC", "Reg-Simple",...
            "Pen-$\ell_{\infty}$-LPEC","Pen-$\ell_{1}$-LPEC",...
            "Feasibility-$\ell_1$", "Feasibility-$\ell_{\infty}$"];
        filename = 'macmpec';
        N_plot = [1,2,4,5,6];
        plot_settings.relative_phase_i = 1;
        plot_settings.absolute_phase_i = 1;
        plot_settings.solved_in_phase_i = 1;
        dtable_reg_simple = dtable(dtable.solver_name == solver_names{2},:);
        max(dtable_reg_simple.n_nlp_total(dtable_reg_simple.success ==1))
        dtable_reg_simple.n_biactive(dtable_reg_simple.n_nlp_total == 15)
        dtable_reg_lpec = dtable(dtable.solver_name == solver_names{1},:);
        max(dtable_reg_lpec.n_nlp_total(dtable_reg_lpec.success ==1))


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

%% Plot results

plot_mpec_benchmark_result(dtable_plot,filename,plot_settings)

%%
% 
% dtable_jj = dtable_plot(dtable_plot.solver_name == solver_names{2},:);
% max(dtable_jj.n_nlp_total(dtable_jj.success == 1))
% hh = dtable_jj.success == 0;
% dtable_jj.problem_name(find(hh));
% % dtable_jj.problem_name(find(hh))
% 
% 
% dtable_ii = dtable_plot(dtable_plot.solver_name == solver_names{4},:);
% hh = dtable_ii.success == 1 & dtable_ii.multiplier_based_stationarity ~= 'S';
% dtable_ii.problem_name(find(hh));
% find(hh')
%% LPEC failure
% dtable_ii = dtable_plot(dtable_plot.solver_name == solver_names{1},:);
% dtable_jj = dtable_plot(dtable_plot.solver_name == solver_names{2},:);
% dtable_jj.problem_name(dtable_ii.success ==1 & dtable_jj.success ==0)


