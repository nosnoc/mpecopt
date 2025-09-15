%% run subseqeuntlly benchmarks of different size
try
    run_nonlinear_mpec_benchmark_small
catch
end
%%
try
    run_nonlinear_mpec_benchmark_medium
catch

end

%% run subseqeuntlly benchmarks of different size
try
    run_nonlinear_mpec_benchmark_large_1
catch
end


%% run subseqeuntlly benchmarks of different size
try
    run_nonlinear_mpec_benchmark_large_2
catch
end