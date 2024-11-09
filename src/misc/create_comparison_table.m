function output_tab = create_comparison_table(dtable,criterion_cell)

output_tab = table;
solver_tab = pivot(dtable, Rows=["solver_name"]);

for ii=1:height(solver_tab)
    for jj=1:length(criterion_cell)
        solver_name = solver_tab.solver_name(ii);
        if ii == 1 && jj == 1
            output_tab.problem_name = char(dtable.problem_name(dtable.solver_name == solver_name));
        end
        data_ii = dtable.(criterion_cell{jj})(dtable.solver_name == solver_name);
        solver_name_str = char(solver_name);
        solver_name_str = regexprep(solver_name_str , ' ', '_');
        solver_name_str = [solver_name_str,'_',criterion_cell{jj}];
        output_tab.(solver_name_str) = data_ii;
    end
end
end