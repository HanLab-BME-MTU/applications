function loadCamera(obj)
if obj.displayEnabled
    % Load the camera properties
    obj.imarisApp.mSurpassCamera.SetOrientationAxisAngle(obj.cameraOrientation(1),obj.cameraOrientation(2),obj.cameraOrientation(3),obj.cameraOrientation(4));
    obj.imarisApp.mSurpassCamera.SetPosition(obj.cameraPosition(1),obj.cameraPosition(2),obj.cameraPosition(3));
    obj.imarisApp.mSurpassCamera.mFocus = obj.cameraFocus;
    obj.imarisApp.mSurpassCamera.mHeight = obj.cameraHeight;
end
end