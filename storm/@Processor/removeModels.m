function removeModels(obj,idx)

idxKeep = setdiff(1:numel(idx),find(idx));

if ~isempty(obj.data.modelType) && ...
        ~isempty(obj.data.modelLength) && ...
        ~isempty(obj.data.modelRes) && ...
        ~isempty(obj.data.modelProj) && ...
        ~isempty(obj.data.modelBezCP) && ...
        ~isempty(obj.data.modelVar) && ...
        ~isempty(obj.data.modelIsOutOfDate)
    
    obj.data.modelType = obj.data.modelType(idxKeep);
    obj.data.modelLength = obj.data.modelLength(idxKeep);
    obj.data.modelRes = obj.data.modelRes(idxKeep);
    obj.data.modelProj = obj.data.modelProj(idxKeep);
    obj.data.modelBezCP = obj.data.modelBezCP(idxKeep);
    obj.data.modelVar = obj.data.modelVar(idxKeep);
    obj.data.modelIsOutOfDate = obj.data.modelIsOutOfDate(idxKeep);
    
end

if ~isempty(obj.data.clusterColor)
    obj.data.clusterColor = obj.data.clusterColor(idxKeep,:);
end

end

