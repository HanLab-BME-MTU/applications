function saveCamera(obj)
if obj.displayEnabled
    % Save the camera properties
    [posX,posY,posZ,angle] = obj.imarisApp.mSurpassCamera.GetOrientationAxisAngle();
    obj.cameraOrientation = [posX,posY,posZ,angle];
    [posX,posY,posZ] = obj.imarisApp.mSurpassCamera.GetPosition();
    obj.cameraPosition = [posX,posY,posZ];
    obj.cameraFocus = obj.imarisApp.mSurpassCamera.mFocus;
    obj.cameraHeight = obj.imarisApp.mSurpassCamera.mHeight;
end
end