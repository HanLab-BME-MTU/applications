classdef Config < handle
    
    % ---------------------
    % This class contains the configuration parameters of the algorithm.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'public',SetAccess = 'public')
        
        % === Config
        configName; % The name of this configuration
        
        % === Debug
        snapshotsEnabled = false; % Enable/disable Imaris snapshots of the current algorithm state (Windows only)
        displayEnabled = false; % Enable/disable the visualization of the intermediate states in Imaris (Windows only)
        snapshotsPath = 'C:\Users\PB93\Desktop\Snapshots\'; % The path where Imaris will store the snapshots
        intermediateResultsTimerEnabled = false; % Enable/disable writing intermediate timer results to the disk
        intermediateResultsDataEnabled = false; % Enable/disable writing intermediate data results to the disk
        
        % === Data
        path; % Location of the raw data
        fileName; % File name of the raw data
        errorX; % Localization precision convolved with model width in X [nm]
        errorY; % in Y [nm]
        errorZ; % in Z [nm]
        roiPosition = [0 0 0]; % Relative ROI position
        roiSize = [0 0 0]; % The size of the ROI 
        
        % === Preprocess
        dataReductionEnabled = false; % Enable/disable the data reduction step
        nReductionRun = 0; % The data reduction factor (Each run reduces the data by about a factor 2)
        reductionEdgeRadius; % Longest edge considered in reduction algorithm
        densityFilteringEnabled = false; % Enable/disable the density based filtering
        neighborBallRadius; % Size of the neighborhood in which the density will be computed
        nNeighborsThreshold; % Threshold on the number of points in the neighborhood
        edgeWidthInitFree = 0; % e <= 0: Disabled, e > 0: Enabled
        subsampleFraction = 1; % 1 == No subsampling, < 1 Fraction of points remaining, > 1 max number of points remaining
        
        % === Orientation detector
        filterLength; % Length of the filer support
        angularSampling = 5; % Step size of the angular sampling
        
        % === Model fitting
        maxDegreeBezier = 3; % Maximal curve complexity
        maxCurvature; % Maximal model curvature (standard fit) or the beta (snakes fit) depending on fitMethod
        betaVar = 10000; % Shape parameter of the variance prior
        modeVar = 1; % Expected mean square of the normalized orthogonal distances
        fitMethod = 1; % 1:Standard 3D, 2:Standard 2D, 3:Snakes 3D, 4:Snakes 2D
        
        % === Edges
        initialEdgeRadiusGeom; % Longest possible edge in the initialization phase
        initialEdgeRadius; % Longest possible edge in the main 
        
        % === Geometric clustering
        nIterGeomMatching; % Number of iterations of the geometric merging
        modelLength; % Model length used to compute the distance criterion
        angleThreshold; % Upper bound on the angle between the orientation of two points allowed to merge
        
        % === Main
        nIterEM; % Maximum expectation-maximization algorithm iterations
        maxIterMerge; % Maximum merging iterations
        
        % === Noise
        nSigmaThreshold; % Multiple of the component standard deviation used to classify points as noise     
                
    end
    
    methods(Static = true)
        obj = load(fullPath);
        obj = read(fullPath);
    end
    
    methods
        function obj = Config()
        end
        
        save(obj,fullPath);
        display(obj);
                
    end
    
end


