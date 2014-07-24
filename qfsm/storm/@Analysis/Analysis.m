classdef Analysis < handle
    
    % ---------------------
    % Class providing the tools needed to analyze the reconstructed model.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'protected',SetAccess = 'protected')
        data;
    end
    
    methods
        function obj = Analysis(data)
            if nargin>0
                obj.data = data;
            end
        end
                
        branches(obj,myImarisVisualizer,edgeRadius,minAngleThreshold,maxAngleThreshold,maxDistanceThreshold);
        bundles(obj,myImarisVisualizer,edgeRadius,distanceThreshold,minOverlapThreshold);
        [meanDist,fragmentation,crosslinking,falsePositivesFraction] = reconstructionScore(obj,sampleInterval,edgeRadius);
        meanSigmaDist = width(obj,display);
        [medianSpacing,probDensFct] = spacings(obj,display);
        maxGap(obj);
        maxCurvature(obj);
        statistics(obj);
        avLen = length(obj);
        density(obj);
        clusterSize(obj);
        meanIntensity = intensity(obj);
        or = orientation(obj,movingAverageSpan,displayText,maxHistValue,zoom)
        [sigX,sigY,sigZ] = residualsXYZ(obj,modelIdx);
        rms = residuals(obj,modelIdx);
        
    end
    
end