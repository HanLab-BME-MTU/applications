%% calculates magnitude, temporal and spatial derivatives
function [magTimeSpace] = whGetDerivatives(features,meanFeats,stdFeats,derivCoeff)
features = features';
[nObs mFeats] = size(features);
curFeatsNorm = (features - repmat(meanFeats,[nObs 1])) ./ repmat(stdFeats,[nObs 1]);

mag = curFeatsNorm * derivCoeff.mag;
timeDeriv = curFeatsNorm * derivCoeff.timeDerive;
spaceDeriv = curFeatsNorm * derivCoeff.spaceDeriv;
magTimeSpace = [mag,timeDeriv,spaceDeriv];
end