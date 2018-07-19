function [minPoint,maxPoint,length,idx,tNormalized] = projectPointsOntoLine(points,startPoint,directionVector)

nPoints = size(points,1);

helpVector = points-repmat(startPoint,nPoints,1);

t = dot(repmat(directionVector,nPoints,1),helpVector,2)/norm(directionVector)^2;

[t,idx] = sort(t);

minPoint = startPoint+t(1)*directionVector;
maxPoint = startPoint+t(end)*directionVector;

length = t(end)-t(1)*norm(directionVector);

tNormalized = (t-t(1))/(t(end)-t(1));

end

