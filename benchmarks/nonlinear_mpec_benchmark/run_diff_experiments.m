% %% run subseqeuntlly benchmarks of different size
% try
%     run_nonlinear_mpec_benchmark_small
% catch
% end
% %%
% try
%     run_nonlinear_mpec_benchmark_medium
% catch
% 
% end
% 
% %% run subseqeuntlly benchmarks of different size
% try
%     run_nonlinear_mpec_benchmark_large_1
% catch
% end


%% run subseqeuntlly benchmarks of different size
try
    run_nonlinear_mpec_benchmark_large_2
catch
end

%% 
try 
    cd ..\macmpec\
    run_macmpec_experiments_general
catch
end