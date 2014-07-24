function weight = geomClusterWeight(obj,p1,p2,meanLocPrec,modelLength,angleThreshold)

% Copy data
point1 = obj.data.points(p1,:); % Inner point
point2 = obj.data.points(p2,:);

orientation1 = obj.data.orientation(p1,:); % Corresponding normalized orientation
orientation2 = obj.data.orientation(p2,:);

% Compute the angle between the segments
if orientation1 == orientation2 % Parallel    
    cosAngle = 1;
else    
    cosAngle = abs(dot(orientation1,orientation2)/(norm(orientation1)*norm(orientation2)));
end

% Compute the distance between the segments
p1 = (point1-modelLength/2*orientation1);
p2 = (point1+modelLength/2*orientation1);
p3 = (point2-modelLength/2*orientation2);
p4 = (point2+modelLength/2*orientation2);

if meanLocPrec(3) == 0 % 2D case
    p1(1:2) = p1(1:2)./meanLocPrec(1:2);
    p2(1:2) = p2(1:2)./meanLocPrec(1:2);
    p3(1:2) = p3(1:2)./meanLocPrec(1:2);
    p4(1:2) = p4(1:2)./meanLocPrec(1:2);
else % 3D case
    p1 = p1./meanLocPrec;
    p2 = p2./meanLocPrec;
    p3 = p3./meanLocPrec;
    p4 = p4./meanLocPrec;
end

% The distance is weighted by meanLocPrec
dist = segments_dist_3d(p1',p2',p3',p4');

% Diable the edges which form a too big angle
if cosAngle < cosd(angleThreshold)
    cosAngle = -1;
end

% Assemble the weight
weight = exp(-0.5*dist.^2)*cosAngle;
   
end