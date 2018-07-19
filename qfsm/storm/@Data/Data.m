classdef Data < handle & matlab.mixin.Copyable
    
    % ---------------------
    % This class contains the input data of the algorithm, intermediate 
    % and final outputs. 
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'public',SetAccess = 'public')
        
        % #####################
        % Data
        % #####################
        
        % ROI
        roiPosition = []; % Absolute position of the region of interest
        roiSize = []; % The size of the region of interest

        % Points
        rawPoints; % The points before preprocessing
        intensity; % Intensity of the point
        frame; % Frame in which the point has been localized
        error; % Localization precision (standard deviation) for x,y and z
        points; % The points after preprocessing
        
        % #####################
        % Outputs
        % #####################
        
        % Points
        orientation; % Normalized vector representing the orientation of a point
        magnitude; % Magnitude of the orientation detection
        neighbors; % Neighbor points indices
        
        % Clusters
        clusters; % Cell array containing point indices
        nullCluster; % Indices of the points belonging to the noise cluster
        clusterColor; % The color of the cluster (Visualization)
        
        % Cluster models
        modelIsOutOfDate; % Boolean: True if the cluster members changed and the model has not yet been updated
        modelType; % The degree of the Bezier curve
        modelLength;
        modelRes; % Residuals of the fit ( n x 3 )
        modelProj; % Projections of the point onto the model (t)
        modelBezCP; % Bezier curve model control points
        modelVar; % The mean square of the normalized distances 
        
        % Edges
        initialEdges; % All allowed edges
        edges; % Current edges
        weights; % Corresponding weights
        
        % Running time
        runTime; % Total algorithm runtime
        
        % #####################
        % Visualization
        % #####################
        
        edgesHistory;
        clustersHistory;
        modelsHistory;
        
        % #####################
        % Simulation
        % #####################
        
        simModelBezCP; % Control points of the model that generated the points
        
    end
    
    properties (GetAccess = 'public',SetAccess = 'private',Dependent = true)
        nPoints;
        nClusters;
        clusterSize;
        parents; % Parent cluster
    end
    
    methods(Static = true)
        obj = load(fullPath);
        obj = read(fullPath);
    end
    
    methods
        % Constructor
        function obj = Data()
        end
        
        display(obj);
        save(obj,fullPath);
        mergeWith(obj,dataToMerge);

        function value = get.nPoints(obj)
            value = size(obj.points,1);
        end
        
        function value = get.nClusters(obj)
            value = size(obj.clusters,1);
        end
        
        function value = get.clusterSize(obj)
            value = cellfun(@numel,obj.clusters);
        end
        
        function parents = get.parents(obj)
            % Update the cluster parents
            nPoints = cellfun(@numel,obj.clusters);
            lastPoint = cumsum(nPoints);
            clusterIdxs = (1:size(lastPoint,1))';
            clustersVector = horzcat(obj.clusters{:})';
            
            clusterIdxs = arrayfun(@(a,b) ones(b,1)*a,clusterIdxs,nPoints,'UniformOutput',false);
            clusterIdxs = vertcat(clusterIdxs{:});

            parents = zeros(obj.nPoints,1);
            parents(clustersVector,1) = clusterIdxs;
        end
        
    end % methods
    
end % class


