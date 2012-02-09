classdef ROISelector < handle
    
    % ---------------------
    % This class enables you to create a preview of the STORM data and it
    % provides a GUI to select a ROI of the data set. The parameters of the
    % ROI are stored as properties.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties
        extentX;
        extentY;
        minX;
        minY;
        previewScaleFactor;
        preview;
        pathName;
        fileName;
        roiPos;
        roiSize;
        
        repetitions = [2,2];
        offset = [1000 1000];
        
    end
    
    methods (Static = true)
        obj = load(fullPath);
    end
    
    methods
        function obj = ROISelector(pathName,fileName)
            if nargin>0
                obj.pathName = pathName;
                obj.fileName = fileName;
            end
        end
        
        run(obj);
        buildPreview(obj);
        h = setupGUI(obj);
        save(obj);
       
    end
    
end


