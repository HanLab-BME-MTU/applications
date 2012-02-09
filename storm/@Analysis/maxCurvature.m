function maxCurvature(obj)
% Compute the maximum curvature of every model
maxCurv = cellfun(@maxCurvatureBezier,obj.data.modelBezCP);

% Display histogram
nBins = round(sqrt(numel(maxCurv)));
hist(maxCurv,nBins);
xlabel('Maximum curvature [nm]');
ylabel('Number of models');

fprintf('Analyzer: Maximum curvature: %f [1/nm]\n',max(maxCurv));
end