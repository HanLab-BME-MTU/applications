function [ind, distance] = distance_two_curves(c1_yx, c2_yx)
%% Calculate the distance of two curve segments( defined by the points on the curve)

% the out put ind is for which points is giving the closest distance
% Make one in columns, one in rows

% Liya Ding

% Default input format
% [ x1 y1
%   x2 y2
%   x3 y3
%   ......
%   xn yn];


if(size(c1_yx,1)==2 && size(c1_yx,2)==2)
    % in the case of 2*2, so assume it us in the needed format    
else    
    if size(c1_yx,1) < size(c1_yx,2)
        c1_yx = c1_yx';
    end
    
    if size(c1_yx,1)==2 && size(c1_yx,2)==1
        c1_yx = c1_yx';
    end
end

if(size(c2_yx,1)==2 && size(c2_yx,2)==2)
    % in the case of 2*2, so assume it us in the assumed format
     c2_yx = c2_yx';
else
    if size(c2_yx,1) > size(c2_yx,2)
        c2_yx = c2_yx';
    end
    
    if size(c2_yx,1)==1 && size(c2_yx,2)==2
        c2_yx = c2_yx';
    end
end

c1_y = c1_yx(:,1);
c1_x = c1_yx(:,2);

c2_y = c2_yx(1,:);
c2_x = c2_yx(2,:);

c1_y_matrix = repmat(c1_y, 1,length(c2_y));
c1_x_matrix = repmat(c1_x, 1,length(c2_y));


c2_y_matrix = repmat(c2_y, length(c1_y),1);
c2_x_matrix = repmat(c2_x, length(c1_y),1);


diff_matrix = (c1_x_matrix - c2_x_matrix).^2 + (c1_y_matrix - c2_y_matrix).^2;

[inda, indb] = find(diff_matrix == min(min(diff_matrix)));

inda = inda(1);
indb = indb(1);

ind = [inda indb];

distance = sqrt(diff_matrix(inda, indb));



