function var = modelDistanceVarianceWeightedWithPrior(res,weights,betaVar,modeVar)

sumWeightedSqDist = sum(sum((res.*weights).^2));
nPoints = size(res,1);

var = (sumWeightedSqDist+2*betaVar)/(2*(betaVar-modeVar)/modeVar+nPoints+2);

end
