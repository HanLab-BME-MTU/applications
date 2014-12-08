
% TODO: use mesh grid
function [cogX, cogY] = cog(bw)
sumX = 0; sumY = 0; num = 0;
for y = 1 : size(bw,1)
    for x = 1 : size(bw,2)
        if bw(y,x) > 0
            sumX = sumX + x;
            sumY = sumY + y;
            num = num + 1;                        
        end
    end
end
cogX = round(sumX/num);
cogY = round(sumY/num);
end