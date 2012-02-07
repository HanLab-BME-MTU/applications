function modelGroundTruth(obj)

nModels = numel(obj.data.simModelBezCP);

if nModels
    pointsStart = zeros(1,3);
    pointsEnd = zeros(1,3);
    
    nCurveSamples = 100;
        
    for m=1:nModels
        modelType = size(obj.data.simModelBezCP{m},1)-1;
        if modelType == 1;
            pointsStart = [pointsStart;obj.data.simModelBezCP{m}(1,:)];
            pointsEnd = [pointsEnd;obj.data.simModelBezCP{m}(2,:)];
        else
            curvePoints = renderBezier(obj.data.simModelBezCP{m},(linspace(0,1,nCurveSamples))');
            pointsStart = [pointsStart;curvePoints(1:(end-1),:)];
            pointsEnd = [pointsEnd;curvePoints(2:end,:)];
        end
    end
    
    pointsStart = pointsStart(2:end,:);
    pointsEnd = pointsEnd(2:end,:);
    obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Ground Truth');
   
else
    disp('Display: The truth is unknown ;-)');
end

end
