function [models res] = getSegmentModels(allFeatures, nFeatures)

nSegments = numel(nFeatures);

models = zeros(nSegments, 4);
res = cell(nSegments, 1);

% pFirst and pLast are indexing allFeatures
pLast = cumsum(nFeatures);
pFirst = pLast-nFeatures+1;

for iSegment = 1:nSegments
    params = allFeatures(pFirst(iSegment):pLast(iSegment), :);
    params = num2cell(params,1);
    [x, y, sx, t] = params{:};
    ct = cos(t);
    st = sin(t);
    x = [x + sx .* ct; x - sx .* ct];
    y = [y + sx .* st; y - sx .* st];
    
    [models(iSegment, :), res{iSegment}] = getSegmentModel(x,y);
end
