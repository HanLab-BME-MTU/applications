function keepWindowOpen(obj)
if obj.displayEnabled
    % If the imaris application handle is deleted in MATLAB, Imaris is kept open
    obj.imarisApp.mUserControl = 1;
end
end
