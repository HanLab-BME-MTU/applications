function logL = logLikelihoodModelAlongAway(res,weights,sigmaDistance,modelLength)

sumWeightedSqDist = sum(sum((res.*weights).^2));
nPoints = size(res,1);

% Curve model
logL = -nPoints*log(sqrt(2*pi)*sigmaDistance*modelLength)-1/(2*sigmaDistance^2)*sumWeightedSqDist;

end
