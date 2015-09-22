function takeSnapshotAndResetScene(obj)
if obj.displayEnabled
    if obj.snapshotsEnabled
        obj.takeSnapshot();
    end
    obj.resetScene();
end
end

