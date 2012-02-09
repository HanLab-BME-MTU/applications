function run(obj)

if exist([obj.pathName obj.fileName(1:end-4) '.prv'],'file')
    data = obj.load([obj.pathName obj.fileName]);
    obj.extentX = data.extentX;
    obj.extentY = data.extentY;
    obj.minX = data.minX;
    obj.minY = data.minY;
    obj.previewScaleFactor = data.previewScaleFactor;
    obj.preview = data.preview;
    disp('ROISelector: Preview read from file!');
else
    obj.buildPreview();
    obj.save();
end

obj.setupGUI();

disp('ROISelector: Done!');

end

