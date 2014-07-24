function buildPreview(obj)

% Read data set
myData = Data.read([obj.pathName obj.fileName]);

disp('ROISelector: Generating preview ...');

% Max size of the preview image
maxPreviewSize = min([round(sqrt(size(myData.points,1))*1.5) 1000]); % 1.5 is userdefined

% Determine size of the preview image
obj.minX = min(myData.points(:,1));
obj.minY = min(myData.points(:,2));
obj.extentX = max(myData.points(:,1))-obj.minX;
obj.extentY = max(myData.points(:,2))-obj.minY;

if obj.extentX < obj.extentY
    sizeY = maxPreviewSize;
    obj.previewScaleFactor = (maxPreviewSize-1)/obj.extentY;
    sizeX = round(obj.previewScaleFactor*obj.extentX)+1;
else
    sizeX = maxPreviewSize;
    obj.previewScaleFactor = (maxPreviewSize-1)/obj.extentX;
    sizeY = round(obj.previewScaleFactor*obj.extentY)+1;
end

obj.preview(1:sizeY,1:sizeX) = 0;

% Translate and scale the data and round to the nearest neighbour
x = round((myData.points(:,1)-obj.minX)*obj.previewScaleFactor)+1;
y = round((myData.points(:,2)-obj.minY)*obj.previewScaleFactor)+1;

% Add the molecules to a 2D array
for i=1:size(x,1)
    obj.preview(y(i),x(i)) = obj.preview(y(i),x(i)) + 1;
end

% Attenuate high intensity points
threshold = quantile(obj.preview(:),.995);
obj.preview(obj.preview > threshold) = threshold;

disp('ROISelector: Preview generated!');

end

