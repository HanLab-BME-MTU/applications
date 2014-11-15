classdef Detections < hgsetget
    % Data encapsulator for detections
    properties (SetAccess = protected)
       xyz; dxyz;    % 3xN or 2xN coordinate
       I; dI;        % Image intensity
       A; dA;        % Image amplitude (when estimated) 
       label;        % An image representing the label of the proposed object.
       MD;           % associated MovieData
       frameIdx;     % Idx in the associated movieData
    end
    methods
        
        function showLabel(vol)
            ip=inputParser;
            ip.add
    


classdef MovieDetections < hgsetget
    % Data encapsulator for detections in a movie
    
    properties (SetAccess = protected)
       detectionsList;      
    end
    
    methods
        %% Constructor
        function obj = MovieDetections(detectionsCell)
           obj.detectionsList=detectionsCell;
        end
    end
    methods (Static,Abstract)
        getPathProperty()
        getFilenameProperty()
    end
end

