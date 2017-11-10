function stop_num = common_str_stop(input_strs)

if length(input_strs)<=1
    stop_num =1;
    return;
end

if iscell(input_strs)==0
    stop_num =1;
    return;
end

reached = 0;

for stop_num = 1 : length(input_strs{1})    
    for i_str = 2 : length(input_strs)
        try
            if(~strcmp(input_strs{1}(1:stop_num),input_strs{i_str}(1:stop_num)))
                reached=1;
                break;
            end
        catch
            reached=1;
            break;
        end
    end
    if(reached==1)
        break;
    end
end