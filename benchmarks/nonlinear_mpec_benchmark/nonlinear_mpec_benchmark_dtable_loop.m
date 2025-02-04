%% Create data struct
dstruct = struct;
dstruct.solver_name = [];
dstruct.problem_name = [];
dstruct.success = logical([]);
dstruct.success_phase_i = logical([]);
dstruct.comp_res = [];
dstruct.cpu_time = [];
dstruct.cpu_time_phase_i = [];
dstruct.cpu_time_phase_ii = [];

dstruct.cpu_time_nlp = [];
dstruct.cpu_time_nlp_phase_i = [];
dstruct.cpu_time_nlp_phase_ii = [];

dstruct.cpu_time_lpec_phase_i = [];
dstruct.cpu_time_lpec_phase_ii = [];
dstruct.cpu_time_lpec = [];
dstruct.n_biactive = [];
dstruct.active_set_changes = [];
dstruct.n_nlp_total = [];
dstruct.n_lpec_total = [];
dstruct.multiplier_based_stationarity = [];
dstruct.b_stationarity = [];
dstruct.solved_in_phase_i = [];
dstruct.max_iterations_reached = [];
dstruct.max_cpu_time_nlp  = [];
dstruct.max_cpu_time_nlp_phase_i  = [];
dstruct.max_cpu_time_nlp_phase_ii  = [];
dstruct.max_cpu_time_lpec = [];
dstruct.problem_infeasible = [];
dstruct.f_lpec = [];

% dstruct.solver_message = [;]
dstruct.prob_num = [];
dstruct.f = [];
n_problems = length(mpecs);
% n_problems = 3;

n_problems = length(mpecs);
n_mpecs = 1:length(mpecs);
n_mpec = 1;
%% Run benchmark
for ii = N_experiments
    name = solver_names(ii);
    options = opts{ii};
    j = 1;
    total_success = 0;
    for jj = n_mpecs
        mpec = mpecs(jj);
        w = mpec.w;
        f = mpec.f_fun(mpec.w);
        g = mpec.g_fun(mpec.w);
        G = mpec.G_fun(mpec.w);
        H = mpec.H_fun(mpec.w);
        w0 = mpec.w0;
        lbw = mpec.lbw;
        ubw = mpec.ubw;
        lbg = mpec.lbg;
        ubg = mpec.ubg;
        mpec_name = mpec.name;
    
        disp([ num2str(jj) '/' num2str(length(mpecs)) ': ' 'solving ' char(mpec.name) ' with solver ' char(name)]);
        fprintf('Problem info:name = %s, n_w = %d, n_g = %d, n_comp = %d\n',name, length(w),length(g),length(G))
        mpec_struct = struct('x',w,'f',f,'g',g,'G',G,'H',H);
        solver_initalization = struct('x0', w0, 'lbx',lbw, 'ubx',ubw,'lbg',lbg,'ubg',ubg);
        % to add non-mpecsol solver you can add ifs here
        [result,stats] = solver_functions{ii}(mpec_struct,solver_initalization,options);

        dstruct.solver_name = [dstruct.solver_name; string(name)];
        dstruct.problem_name = [dstruct.problem_name; string(mpec_name)];
        dstruct.success = [dstruct.success; logical(stats.success)];
        dstruct.success_phase_i = [dstruct.success_phase_i; logical(stats.success_phase_i)];

        dstruct.comp_res = [dstruct.comp_res; stats.comp_res];
        dstruct.n_nlp_total = [dstruct.n_nlp_total; stats.n_nlp_total];
        dstruct.n_lpec_total = [dstruct.n_lpec_total; stats.n_lpec_total];

        dstruct.max_cpu_time_nlp = [dstruct.max_cpu_time_nlp; max([stats.iter.cpu_time_nlp_phase_i_iter(:);stats.iter.cpu_time_nlp_phase_ii_iter(:)])];
        dstruct.max_cpu_time_nlp_phase_i = [dstruct.max_cpu_time_nlp_phase_i; max(stats.iter.cpu_time_nlp_phase_i_iter)];
        dstruct.max_cpu_time_nlp_phase_ii = [dstruct.max_cpu_time_nlp_phase_ii; max(stats.iter.cpu_time_nlp_phase_ii_iter)];
        dstruct.cpu_time = [dstruct.cpu_time; stats.cpu_time_total];
        dstruct.cpu_time_phase_i = [dstruct.cpu_time_phase_i; stats.cpu_time_phase_i];
        dstruct.cpu_time_phase_ii = [dstruct.cpu_time_phase_ii; stats.cpu_time_phase_ii];

        dstruct.cpu_time_nlp = [dstruct.cpu_time_nlp; stats.cpu_time_nlp];
        dstruct.cpu_time_nlp_phase_i = [dstruct.cpu_time_nlp_phase_i; stats.cpu_time_nlp_phase_i];
        dstruct.cpu_time_nlp_phase_ii = [dstruct.cpu_time_nlp_phase_ii; stats.cpu_time_nlp_phase_ii];

        dstruct.max_cpu_time_lpec = [dstruct.max_cpu_time_lpec; max([stats.iter.cpu_time_lpec_phase_i_iter';stats.iter.cpu_time_lpec_phase_ii_iter'])];
        dstruct.cpu_time_lpec = [dstruct.cpu_time_lpec; stats.cpu_time_lpec];
        dstruct.cpu_time_lpec_phase_i = [dstruct.cpu_time_lpec_phase_i; stats.cpu_time_lpec_phase_i];
        dstruct.cpu_time_lpec_phase_ii = [dstruct.cpu_time_lpec_phase_ii; stats.cpu_time_lpec_phase_ii];       
        
        dstruct.n_biactive = [dstruct.n_biactive; stats.n_biactive];
        dstruct.f_lpec = [dstruct.f_lpec; stats.f_lpec];

        dstruct.active_set_changes = [dstruct.active_set_changes; sum(stats.iter.active_set_changes)];
        % dstruct.success_despite_infeasiblity = [dstruct.success_despite_infeasiblity ;stats.success_despite_infeasiblity];

        dstruct.solved_in_phase_i = [dstruct.solved_in_phase_i ;stats.solved_in_phase_i];
        dstruct.max_iterations_reached = [dstruct.max_iterations_reached ; stats.max_iterations_reached];
        dstruct.problem_infeasible = [dstruct.problem_infeasible; stats.problem_infeasible];
        dstruct.prob_num = [dstruct.prob_num; jj];
        dstruct.f = [dstruct.f; result.f];
        dstruct.b_stationarity = [dstruct.b_stationarity; stats.b_stationarity];
        dstruct.multiplier_based_stationarity  = [dstruct.multiplier_based_stationarity; stats.multiplier_based_stationarity];
        % j = j+1;
        total_success = total_success+stats.success;
        fprintf(['success rate so far (%d of %d) or %2.2f precent \n\n'],total_success,jj,100*total_success/jj);
        jj = jj+1;

    end
    % save intermediate reuslts 
    dtable1 = struct2table(dstruct);
    save([results_name '_' num2str(ii)],"dtable1");
    % pause(90); % cool down cpu pause
end
%% Check results and plot
dtable = struct2table(dstruct);
save(results_name,"dtable");
