function fitAndSaveCamera(obj)
if obj.displayEnabled
    % Fit the camera view to the scene content
    obj.fitCamera();
    
    % Save the camera properties
    obj.saveCamera();
end
end