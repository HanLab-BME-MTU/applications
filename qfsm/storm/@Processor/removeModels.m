function removeModels(obj,idx)

idxKeep = setdiff(1:numel(idx),find(idx));

if ~isempty(obj.data.modelType)
    obj.data.modelType = obj.data.modelType(idxKeep);
end

if ~isempty(obj.data.modelLength)
    obj.data.modelLength = obj.data.modelLength(idxKeep);
end

if ~isempty(obj.data.modelRes)
    obj.data.modelRes = obj.data.modelRes(idxKeep);
end

if ~isempty(obj.data.modelProj)
    obj.data.modelProj = obj.data.modelProj(idxKeep);
end

if ~isempty(obj.data.modelBezCP)
    obj.data.modelBezCP = obj.data.modelBezCP(idxKeep);
end

if ~isempty(obj.data.modelVar)
    obj.data.modelVar = obj.data.modelVar(idxKeep);
end

if ~isempty(obj.data.modelIsOutOfDate)
    obj.data.modelIsOutOfDate = obj.data.modelIsOutOfDate(idxKeep);
end

if ~isempty(obj.data.clusterColor)
    obj.data.clusterColor = obj.data.clusterColor(idxKeep,:);
end

end

