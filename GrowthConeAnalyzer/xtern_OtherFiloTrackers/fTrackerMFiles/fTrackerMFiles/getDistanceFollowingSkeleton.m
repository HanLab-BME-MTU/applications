function Distance = getDistanceFollowingSkeleton(PositionVectices,shotestPath)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

Distance = 0;
for it = 1:length(shotestPath)-1
   Distance = Distance +  sqrt(sum((PositionVectices(shotestPath(it),:)-PositionVectices(shotestPath(it+1),:)).^2));
end