function saveModelsToHistory(obj)
obj.data.modelsHistory{numel(obj.data.modelsHistory)+1} = obj.data.modelsBezCP;
end