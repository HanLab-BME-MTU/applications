function iCorrespond = correspondingIndices(pointsIn,curve2,maxDist)

if nargin < 3 || isempty(maxDist)
    maxDist = Inf;
end

%Finds the indices for the points in curve 2 which are closest to the input
%points

nPoints = size(pointsIn,2);


for j = 1:nPoints
   
    %Calculate distance to each point on curve 2
    currDists = sqrt( (curve2(1,:) - pointsIn(1,j)).^2 + (curve2(2,:) - pointsIn(2,j)).^2 );
    
    %Find closes point
    [minDist,iTmp] = min( currDists );
    if minDist <= maxDist
        iCorrespond(j) = iTmp;
    else
        iCorrespond(j) = NaN;
    end
    
end