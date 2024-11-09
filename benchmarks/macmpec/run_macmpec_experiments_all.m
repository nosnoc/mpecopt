%% run all experiments
% pause(300)

% try
% clear all;
% % run_macmpec_experiments_phase_i;
% catch
% end

pause(900)
try
clear all;
run_macmpec_experiments_lpec;
catch
end

try
run_macmpec_experiments_general;
catch
end

pause(1800)
try
cd ..\nonlinear_mpec_benchmark\
run_nonlinear_mpec_benchmark
catch
end

