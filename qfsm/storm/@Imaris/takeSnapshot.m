function takeSnapshot(obj)
if obj.displayEnabled
    if obj.snapshotsEnabled
        disp('Imaris: Cheeeeeeeese!');
        obj.snapshotsCounter = obj.snapshotsCounter+1;
        obj.loadCamera();
        obj.imarisApp.SaveSnapShot([obj.snapshotsPath num2str(obj.snapshotsCounter) '.tif']);
    end
end
end