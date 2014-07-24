function lineWithDensityGradient(obj,startPoint,endPoint,densityStart,densityEnd)
orientation = endPoint-startPoint;
length = norm(orientation);

if strcmp(obj.samplingType,'random')
    meanDensity = (densityStart+densityEnd)/2;
    nSamples = round(meanDensity*length);
    
    % Compute probability function
    n = 100*nSamples;
    x = linspace(densityEnd,densityStart,n);
    m = norm(endPoint-startPoint)/(densityEnd-densityStart);
    q = norm(startPoint)-m*densityStart;
    p = m*x + q;
    p = p/sum(p);
    
    % Sample the probability function
    uni = rand(nSamples,1);
    cumprob = [0 cumsum(p)];
    samplePos = zeros(nSamples,1);
    for j=1:n
        ind = (uni>cumprob(j)) & (uni<=cumprob(j+1));
        samplePos(ind) = (j-1)/(n-1);
    end
elseif strcmp(obj.samplingType,'regular')
    samplePos(1) = 0;
    i = 1;
    while samplePos(i) < length
        density = densityStart+(densityEnd-densityStart)*samplePos(i)/length;
        samplePos(i+1,1) = samplePos(i) + 1/density;
        i = i + 1;
    end
    nSamples = numel(samplePos);
    samplePos = samplePos/length;
end

points = repmat(startPoint,nSamples,1)+repmat(orientation,nSamples,1).*repmat(samplePos,1,3);

obj.data.points = [obj.data.points;points];
end