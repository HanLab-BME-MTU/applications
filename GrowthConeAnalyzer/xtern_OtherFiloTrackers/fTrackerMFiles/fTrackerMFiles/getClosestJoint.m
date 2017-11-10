function jointIds=getClosestJoint(xNear, yNear, allJoints)

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


for itNear=1:length(xNear)
    [value, jointIds(itNear)]= min(sqrt(sum((allJoints-[xNear(itNear)*ones(size(allJoints,1),1),yNear(itNear)*ones(size(allJoints,1),1)]).^2,2)));
end

