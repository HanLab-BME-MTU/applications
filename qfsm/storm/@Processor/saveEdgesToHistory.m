function saveEdgesToHistory(obj)
obj.data.edgesHistory{numel(obj.data.edgesHistory)+1} = obj.data.edges;
end