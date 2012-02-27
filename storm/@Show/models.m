function models(obj)

nCurveSamples = 100;

linCurve = obj.data.modelBezCP(obj.data.modelType == 1);
linCurvePointsStart = cellfun(@(a) a(1,:),linCurve,'UniformOutput',false);
linCurvePointsEnd = cellfun(@(a) a(2,:),linCurve,'UniformOutput',false);
linCurvePointsStart = vertcat(linCurvePointsStart{:});
linCurvePointsEnd = vertcat(linCurvePointsEnd{:});

curvePoints = cellfun(@(a) renderBezier(a,(linspace(0,1,nCurveSamples))'),obj.data.modelBezCP(obj.data.modelType >= 2),'UniformOutput',false);
curvePointsStart = cellfun(@(a) a(1:end-2,:),curvePoints,'UniformOutput',false);
curvePointsEnd = cellfun(@(a) a(3:end,:),curvePoints,'UniformOutput',false);
curvePointsStart = vertcat(curvePointsStart{:});
curvePointsEnd = vertcat(curvePointsEnd{:});

curvePointsStart = [curvePointsStart;linCurvePointsStart];
curvePointsEnd = [curvePointsEnd;linCurvePointsEnd];

if ~isempty(curvePointsStart)
    obj.imaris.displaySegments(curvePointsStart,curvePointsEnd,'Display: Models');
end

end