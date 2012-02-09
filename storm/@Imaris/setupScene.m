function setupScene(obj)
if obj.displayEnabled
    % Add a light source
    lightSource = obj.imarisApp.mFactory.CreateLightSource;
    lightSource.mName = 'Imaris: Light Source';
    obj.imarisApp.mSurpassScene.AddChild(lightSource);
    obj.imarisApp.mDataSet.mUnit = 'nm';
end
end
