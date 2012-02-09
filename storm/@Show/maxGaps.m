function maxGaps(obj)

pointsStart = zeros(1,3);
pointsEnd = zeros(1,3);

nCurveSamples = 20;

for c=1:obj.data.nClusters
    
    % Determine the maximum gap of the curve
    t = sort(obj.data.modelProj{c});
    startTs = t(2:end);
    endTs = t(1:end-1);
    fun = @(startTs,endTs) lengthBezier(obj.data.modelBezCP{c},startTs,endTs);
    gaps = arrayfun(fun,startTs,endTs);
    [~,idx] = max(gaps);
    startT = startTs(idx);
    endT = endTs(idx);
    
    % Compute curve points
    maxGapT = linspace(startT,endT,nCurveSamples)';
    curvePoints = renderBezier(obj.data.modelBezCP{c},maxGapT);
    
    pointsStart = [pointsStart;curvePoints(1:(end-1),:)];
    pointsEnd = [pointsEnd;curvePoints(2:end,:)];
    
end

pointsStart = pointsStart(2:end,:);
pointsEnd = pointsEnd(2:end,:);
obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Max Gaps');

end