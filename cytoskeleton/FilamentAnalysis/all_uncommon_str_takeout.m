function output_strs = all_uncommon_str_takeout(input_strs)

% old fashion, when only one file, only the last character is kept.
if iscell(input_strs)==0
    output_strs = input_strs(end);
    return;
end


stop_num = common_str_stop(input_strs);
cutback_num = common_str_cutback(input_strs);



for i_str = 1 : length(input_strs)
    output_strs{i_str} = input_strs{i_str}(stop_num:end-cutback_num);
end
