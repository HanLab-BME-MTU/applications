function fitCamera(obj)
if obj.displayEnabled
    % Fit the camera view to the scene content
    obj.imarisApp.mSurpassCamera.Fit();
end
end