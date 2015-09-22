classdef ForegroundSegmentationModule < handle
% Abstract base class for all foreground segmentation modules

    properties (SetAccess = protected)        
        parameters
    end
     
    % All derived classes must implement these methods
    methods(Abstract)

        % returns a nice descriptive name of the segmentation algorithm
        [strName] = getName(obj);
        
        % runs the algorithm on the imageData provided and returns a
        % binary foreground mask
        [imForegroundMask] = segmentImage(obj, imageData, spacing);
        
        % fires a GUI and allows the user to specify the parameters
        % of the algorithm
        uiSetParameters(obj);

    end

    methods
        
        % set parameters
        function setParameters(obj, parameters)
            
            if ~isempty(setxor(fieldnames(obj.parameters), fieldnames(parameters)))
                error('Invalid parameter sturcture');
            end
            
            obj.parameters = parameters;
            
        end
        
    end
    
end