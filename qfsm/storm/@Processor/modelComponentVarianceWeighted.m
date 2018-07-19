function var = modelComponentVarianceWeighted(res,weights,varMin)

sumWeightedSqDist = sum(sum((res.*weights).^2));
nPoints = size(res,1);

var = 1/(2*nPoints)*sumWeightedSqDist;

var = max(var,varMin);

end

