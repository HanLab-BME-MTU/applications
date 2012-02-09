function line(obj,startPoint,endPoint,nSamples)
if strcmp(obj.samplingType,'random')
    % Random sampling
    points = repmat(startPoint,nSamples,1)+repmat(rand(nSamples,1),1,3).*repmat(endPoint-startPoint,nSamples,1);
elseif strcmp(obj.samplingType,'regular')
    % Regular sampling
    length = endPoint-startPoint;
    sampleStep = length/(nSamples-1);
    points = repmat(startPoint,nSamples,1)+repmat((0:nSamples-1)',1,3).*repmat(sampleStep,nSamples,1);
end
obj.data.points = [obj.data.points;points];
end