function [ aligned, events ] = alignExtrema( extrema , period, unwrap, truncate)
%alignExtrema Align extrema tracks as K decreases (time increases)

if(nargin < 2 || isempty(period))
    period = 2*pi;
end
if(nargin < 3)
    unwrap = false;
end
if(nargin < 4)
    truncate = true;
end

aligned = sort(extrema);
nExtrema = sum(~isnan(aligned));
totalExtrema = max(nExtrema);
aligned = aligned(1:totalExtrema,:);

currentCost = nansum(abs(diff(aligned,1,2)));

% Detect change in number of extrema and crossings of the periodic boundary
events = find(diff(nExtrema) ~= 0 | currentCost > period - currentCost);
% TODO: trigger alignment event if total difference is large due to extrema
% heading over periodic boundary
% events = 1:length(nExtrema)-1;

for e = events
    cost = abs(bsxfun(@minus,aligned(:,e).',aligned(:,e+1)));
%     wrap = cost > period;
%     cost(wrap) = mod(cost(wrap),period);
    cost = min(abs(period-cost),cost);
    max_cost = max(cost(:));
    if(isnan(max_cost))
        % Cost matrix is all NaN
    else
        cost(isnan(cost)) = max(cost(:))+1;
        [link12,link21] = lap(cost);       
        aligned(:,e+1:end) = aligned(link21,e+1:end);
    end
end

if(unwrap)
    aligned = orientationSpace.diffusion.unwrapExtrema(aligned, events, period);
end
if(~truncate)
    extrema(1:totalExtrema,:) = aligned;
    aligned = extrema;
end

end

