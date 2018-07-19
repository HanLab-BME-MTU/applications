function clearScene(obj)
if obj.displayEnabled
    % Get the number of children of the scene
    nChildren = obj.imarisApp.mSurpassScene.GetNumberOfChildren;
    childrenIndexes = flipdim(0:nChildren-1,2);
    
    % Remove all the children
    for child=childrenIndexes
        dataItem = obj.imarisApp.mSurpassScene.GetChild(child);
        obj.imarisApp.mSurpassScene.RemoveChild(dataItem);
    end
end
end