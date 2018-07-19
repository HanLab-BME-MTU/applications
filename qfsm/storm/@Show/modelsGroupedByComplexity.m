function modelsGroupedByComplexity(obj)

nCurveSamples = 100;

linCurve = obj.data.modelBezCP(obj.data.modelType == 1);
linCurvePointsStart = cellfun(@(a) a(1,:),linCurve,'UniformOutput',false);
linCurvePointsEnd = cellfun(@(a) a(2,:),linCurve,'UniformOutput',false);
linCurvePointsStart = vertcat(linCurvePointsStart{:});
linCurvePointsEnd = vertcat(linCurvePointsEnd{:});

quadCurvePoints = cellfun(@(a) renderBezier(a,(linspace(0,1,nCurveSamples))'),obj.data.modelBezCP(obj.data.modelType == 2),'UniformOutput',false);
quadCurvePointsStart = cellfun(@(a) a(1:end-1,:),quadCurvePoints,'UniformOutput',false);
quadCurvePointsEnd = cellfun(@(a) a(2:end,:),quadCurvePoints,'UniformOutput',false);
quadCurvePointsStart = vertcat(quadCurvePointsStart{:});
quadCurvePointsEnd = vertcat(quadCurvePointsEnd{:});

cubCurvePoints = cellfun(@(a) renderBezier(a,(linspace(0,1,nCurveSamples))'),obj.data.modelBezCP(obj.data.modelType == 3),'UniformOutput',false);
cubCurvePointsStart = cellfun(@(a) a(1:end-1,:),cubCurvePoints,'UniformOutput',false);
cubCurvePointsEnd = cellfun(@(a) a(2:end,:),cubCurvePoints,'UniformOutput',false);
cubCurvePointsStart = vertcat(cubCurvePointsStart{:});
cubCurvePointsEnd = vertcat(cubCurvePointsEnd{:});

if ~isempty(linCurvePointsStart)
    obj.imaris.displaySegments(linCurvePointsStart,linCurvePointsEnd,'Display: Linear Models');
end
if ~isempty(quadCurvePointsStart)
    obj.imaris.displaySegments(quadCurvePointsStart,quadCurvePointsEnd,'Display: Quadratic Models');
end
if ~isempty(cubCurvePointsStart)
    obj.imaris.displaySegments(cubCurvePointsStart,cubCurvePointsEnd,'Display: Cubic Models');
end

end