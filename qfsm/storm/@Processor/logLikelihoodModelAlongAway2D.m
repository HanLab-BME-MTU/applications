function logL = logLikelihoodModelAlongAway2D(res,weights,sigmaComponent,modelLength)

sumWeightedSqDist = sum(sum((res.*weights).^2));
nPoints = size(res,1);

if modelLength == 0 % Point model
    logL = -nPoints*3*log(sqrt(2*pi)*sigmaComponent)-1/(2*sigmaComponent^2)*sumWeightedSqDist;
else % Curve model
    logL = -nPoints*log(2*pi*sigmaComponent^2*modelLength)-1/(2*sigmaComponent^2)*sumWeightedSqDist;
end

end

