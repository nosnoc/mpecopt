function [varargout] = unfold_struct(input,mode)
% mode = 'caller'; % copy fields locally in the function where unfold_struct was called
% mode = 'base'; % copy fields  into the main workspace

%% 
import casadi.*
names = fieldnames(input);
for i=1:length(names)
    eval([names{i} '=input.' names{i} ';']);
%     if eval(['~exist( '''  names{i} ''', ''var'')'])
        %         Check if a variable exists in the workspace, within a function
        assignin(mode, names{i}, eval(names{i}));
%     end
end
end




% for i=1:length(names)
%     eval([names{i} '=input.' names{i} ]);
%     eval(['varargout{' num2str(i) '} =' names{i}])
% end
