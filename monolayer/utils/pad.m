function [ipad] = pad(i,digits)

if nargin < 2 || digits < 2
    ipad = sprintf('%d',i);
else
    
    if digits == 3
        i = double(i);
        if (floor(i / 100) > 0)
            ipad = sprintf('%d',i);
        else if (floor(i / 10) > 0)
                ipad = sprintf('0%d',i);
            else
                ipad = sprintf('00%d',i);
            end
        end
    else if digits == 2
            if (floor(i / 10) > 0)
                ipad = sprintf('%d',i);
            else
                ipad = sprintf('0%d',i);
            end
        end
    end
end
end