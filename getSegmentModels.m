function [models res] = getSegmentModels(allFeatures)

nModels = numel(allFeatures);

models = zeros(nModels, 4);
res = cell(nModels, 1);

for iModel = 1:nModels
    params = allFeatures{iModel};
    params = num2cell(params,1);
    [x, y, sx, t] = params{:};
    ct = cos(t);
    st = sin(t);
    x = [x + sx .* ct; x - sx .* ct];
    y = [y + sx .* st; y - sx .* st];
    
    [models(iModel, :), res{iModel}] = getSegmentModel(x,y);
end
