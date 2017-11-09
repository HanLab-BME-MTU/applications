function sample_matrix = boxplot_combine_multi_array(sample_cell)
% combine sample sets with different number of samples, for boxplot

no_sets = numel(sample_cell);

sample_count_array = zeros(no_sets,1);

for i = 1 : no_sets
    sample_count_array(i) = numel(sample_cell{i}(:));
end

sample_matrix = nan(max(sample_count_array),no_sets);

for i = 1 : no_sets
    sample_matrix(1:sample_count_array(i),i) = sample_cell{i}(:);    
end

