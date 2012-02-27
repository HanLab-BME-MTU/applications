classdef Show < handle
    
    % ---------------------
    % The content of a Data object can be visualized with this class
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'private',SetAccess = 'private')
        data; % Data object handle
    end
    
    properties (GetAccess = 'public',SetAccess = 'private')
        imaris; % Imaris handle
    end
    
    properties(Constant = true)
        nullClusterColor = [0.5 0.5 0.5 0]; % Noise cluster color
        pointColor = [0.0 1.0 0.0 0.0];
        projectionColor = [1.0 0.0 0.0 0.0];
        pointSize = 4*0.6*5;% = 4*0.6; % = 4; % The diameter of the points
        projectionSize = 3;
        pointModelSize = 15;
        pointModelColor = [1.0 1.0 1.0 0.9];
        controlPointSize = 15;
        controlPointColor = [0.0 0.0 1.0 0.0];
        nColorSlices = 20;
    end
    
    methods
        function obj = Show(data,imaris)
            if nargin==2
                obj.data = data;
                obj.imaris = imaris;
            elseif nargin==1
                obj.data = data;
                obj.imaris = Imaris();
                obj.imaris.setupScene();
            end
        end
        
        points(obj,varargin);
        pointsRaw(obj);
        pointsWithZColorCoded(obj);
        pointsWithWidth(obj);
        clusters(obj,idx);
        clustersAsChains(obj);
        edges(obj);
        initialEdges(obj);
        nullCluster(obj);
        maxGaps(obj);
        projections(obj);
        residuals(obj);
        models(obj);
        modelsControlPoints(obj);
        modelsGroupedByComplexity(obj);
        orientation(obj,modelLength);
        orientationWithVariableLength(obj,maxModelLength);
        modelGroundTruth(obj);
        clustersHistoryLast(obj);
        edgesHistoryLast(obj,idx);
        roi(obj);
        
    end
    
end

