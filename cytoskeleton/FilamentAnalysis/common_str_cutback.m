function stop_num = common_str_cutback(input_strs)

 stop_num =0;
 
if length(input_strs)<=1
    stop_num =0;
    return;
end

if iscell(input_strs)==0
    stop_num =0;
    return;
end


reached = 0;

for stop_num = 1 : length(input_strs{1})    
    for i_str = 2 : length(input_strs)
        try
            if(~strcmp(input_strs{1}(end-stop_num:end),input_strs{i_str}(end-stop_num:end)))
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