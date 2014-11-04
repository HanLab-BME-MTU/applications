function [line_x, line_y] = line_following_with_limit(bw_in, radius, end_point_x,end_point_y)
% following a line giving the skeleton bw and the starting point (which should be a end point)
% Input: bw_in: the bw image with the line pixels as 1
%        radius: the pixel length of the line wanted; if not giving, then full length
%        end_point_x/y: the end point to start following
% output: line_x/y: the traced line pixels (in order)


% Everything + 1 to avoid trouble of boundary points

if(isempty(radius))
    radius = sum(sum(bw_in))-1;
end

line_x = end_point_x+1;
line_y = end_point_y+1;

bw_working = zeros(size(bw_in,1)+2,size(bw_in,2)+2);
bw_working(2:end-1,2:end-1) = bw_in;

for i = 1 : radius
    if(bw_working(line_y(end),line_x(end))==0)
        % if this is already 0, means there is no more points to following
        % in this case(the loop will be null), so stop
        break;
    end
    
    bw_working(line_y(end),line_x(end))=0;

    current_x = line_x(end);
    current_y = line_y(end);
    
    if(bw_working(current_y-1,current_x)==1)
        
        line_y = [line_y current_y-1];
        line_x = [line_x current_x];
        
        continue;
    end
    
    if(bw_working(current_y,current_x+1)==1)
        
        line_y = [line_y current_y];
        line_x = [line_x current_x+1];
        
        continue;
    end
    
    if(bw_working(current_y+1,current_x)==1)
        
        line_y = [line_y current_y+1];
        line_x = [line_x current_x];
        
        continue;
    end
    
    if(bw_working(current_y,current_x-1)==1)
        
        line_y = [line_y current_y];
        line_x = [line_x current_x-1];
        
        continue;
    end
    
    
    if(bw_working(current_y-1,current_x-1)==1)
        
        line_y = [line_y current_y-1];
        line_x = [line_x current_x-1];

        continue;
    end

    if(bw_working(current_y-1,current_x+1)==1)
        
        line_y = [line_y current_y-1];
        line_x = [line_x current_x+1];

        continue;
    end

    if(bw_working(current_y+1,current_x+1)==1)
        
        line_y = [line_y current_y+1];
        line_x = [line_x current_x+1];

        continue;
    end

    if(bw_working(current_y+1,current_x-1)==1)
        
        line_y = [line_y current_y+1];
        line_x = [line_x current_x-1];

        continue;
    end
        
end

% get rid of the + 1
line_x=line_x-1;
line_y=line_y-1;



            