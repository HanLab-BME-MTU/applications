function [ aligned ] = alignExtrema( extrema )
%alignExtrema Align extrema tracks as K decreases (time increases)

aligned = sort(extrema);
nExtrema = sum(~isnan(aligned));
totalExtrema = max(nExtrema);
aligned = aligned(1:totalExtrema,:);

events = find(diff(nExtrema) ~= 0);

period = 2*pi;

for e = events
    cost = abs(bsxfun(@minus,aligned(:,e).',aligned(:,e+1)));
    wrap = cost > period;
    cost(wrap) = mod(cost(wrap),period);
    cost(isnan(cost)) = max(cost(:));
    [link12,link21] = lap(cost);
    aligned(:,e+1:end) = aligned(link12,e+1:end);
end

end

