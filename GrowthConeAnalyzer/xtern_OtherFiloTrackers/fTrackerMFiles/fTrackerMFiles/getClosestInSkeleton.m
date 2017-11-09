function [xSkel, ySkel]=getClosestInSkeleton(xNear, yNear, skeleton)

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca

    
[row, col]=find(skeleton);
closeDistance=size(skeleton, 1)+size(skeleton, 2);
for it=1:size(row, 1)
    thisDistance=sqrt((row(it)-yNear)^2+(col(it)-xNear)^2);
    if thisDistance<closeDistance
        closeDistance=thisDistance;
        xSkel=col(it);
        ySkel=row(it);
    end
end


        
        