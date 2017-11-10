function [ I ] = drawRadialSegments( angles , startRadius, stopRadius, sz, mode)
%drawRadialSegments draw lines at the angles going symmetrically through the
%center

if(nargin < 4)
    sz = 101;
end
if(nargin < 5)
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

if(nargin < 3)
    stopRadius = floor(I_size(1)/2);
end

% radius = floor(I_size(1)/2);
center = I_size / 2;
% I(center(1),1:center) = 1;
for angle = angles
    unitVector = [cos(angle) sin(angle)];
    startProjection = unitVector*startRadius;
    stopProjection = unitVector*stopRadius;
    p = bresenham(round(center + stopProjection), round(center+startProjection),8);
    idx = sub2ind(I_size, p(:,2), p(:,1));
    switch(mode)
        case 1 % OR
            I(idx) = 1;
        case 2
            I(idx) = I(idx) +1;
    end
end


end

