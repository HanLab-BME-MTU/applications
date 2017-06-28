function [ aligned ] = alignExtrema( extrema , period)
%alignExtrema Align extrema tracks as K decreases (time increases)

if(nargin < 2)
    period = 2*pi;
end

aligned = sort(extrema);
nExtrema = sum(~isnan(aligned));
totalExtrema = max(nExtrema);
aligned = aligned(1:totalExtrema,:);

events = find(diff(nExtrema) ~= 0);
% TODO: trigger alignment event if total difference is large due to extrema
% heading over periodic boundary
% events = 1:length(nExtrema)-1;

for e = events
    cost = abs(bsxfun(@minus,aligned(:,e).',aligned(:,e+1)));
%     wrap = cost > period;
%     cost(wrap) = mod(cost(wrap),period);
    cost = min(period-cost,cost);
    cost(isnan(cost)) = max(cost(:))+1;
    [link12,link21] = lap(cost);
    aligned(:,e+1:end) = aligned(link21,e+1:end);
end

end

