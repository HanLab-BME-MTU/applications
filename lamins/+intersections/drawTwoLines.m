function [ I ] = drawTwoLines( angles, sz, mode)
%drawTwoLines draw lines at the angles going symmetrically through the
%center

if(nargin < 2)
    sz = 101;
end
if(nargin < 3)
    mode = 1;
else
    switch(mode)
        case 'or'
            mode = 1;
        case 'add'
            mode = 2;
        otherwise
            if(~isnumeric(mode))
                warning('intersections.drawRadialLines:Did not understand mode parmeter');
                mode = 1;
            end
    end
end

I = zeros(sz);
I_size = size(I);
radius = I_size(1)/2 -1;
center = I_size / 2;
% I(center(1),:) = 1;
for angle = angles
    projection = [cos(angle) sin(angle)]*radius;
    p = bresenham(round(center + projection), round(center - projection),8);
    idx = sub2ind(I_size, p(:,2), p(:,1));
    switch(mode)
        case 1 % OR
            I(idx) = 1;
        case 2
            I(idx) = I(idx) +1;
    end
end

end

