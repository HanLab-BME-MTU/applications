function [distancia, destination]=getDistance(origin,posibleDestinations)

% Returns the shortest of the distance between the origin and all possible
% destinations

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


for it=1:size(posibleDestinations, 1)
    distances(it) = sqrt((origin(1,1)-posibleDestinations(it, 1))^2+(origin(1,2)-posibleDestinations(it, 2))^2);
end
if exist('distances','var') == 1
    distancia=min(distances);
else
    ''
end
destination=posibleDestinations(find(distances==distancia),:);
destination=[destination(1, 1) destination(1, 2)];