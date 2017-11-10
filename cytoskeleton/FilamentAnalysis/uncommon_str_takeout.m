function output_strs = uncommon_str_takeout(input_strs)

if iscell(input_strs)==0
    output_strs = input_strs;
    return;
end


stop_num = common_str_stop(input_strs);
cutback_num = common_str_cutback(input_strs);


%% stop_num=1 to keep all the head paddings so that less confusion
%  in the ordering of the image seqeucne.

stop_num=1;

for i_str = 1 : length(input_strs)
    output_strs{i_str} = input_strs{i_str}(stop_num:end-cutback_num);
end
