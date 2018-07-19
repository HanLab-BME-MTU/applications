classdef Processor < handle
    
    % ---------------------
    % This class provides all the needed methods to process the data. Its
    % methods are the building blocks of the STORM reconstruction algorithm.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties
        data;
    end
    
    properties (GetAccess = 'private',SetAccess = 'private')
        edgeModelType;
        edgeModelLength;
        edgeModelRes;
        edgeModelProj;
        edgeModelBezCP;
        edgeModelVar;
    end
    
    methods
        function obj = Processor(data)
            if nargin>0
                obj.data = data;
            end
        end
        
        setErrorArray(obj,errorX,errorY,errorZ);
        cropData(obj,roiPosRel,roiSize);
        centerData(obj);
        centerDataRevert(obj,cfg);
        subsamplePoints(obj,limit);
        densityFilter(obj,nNeighborsThreshold,ballRadius);
        nearestNeighborClutterRemoval(obj,k,mode);
        dataReduction(obj,initialEdgeRadius,reductionRuns);
        computeOrientation(obj,filterLength,filterAngularSampling,nBins,minBinResponse);
        [weights,hausDorffEdges] = updateGeomWeights(obj,sigmaFilament,modelLength,angleThreshold);
        mergeClustersGeome(obj,sigmaFilament,modelLength,angleThreshold);
        edgeLength = initEdges(obj,initialEdgeRadius);
        updateEdges(obj);
        updateEdgesEndPoints(obj,dMax);
        updateEdgesAnisotropic(obj,dRef,alpha,samplePeriod,dMaxAlong,dMinAway);
        initClusters(obj);
        initModels(obj,betaVar,modeVar);
        updateModels(obj,maxCurvature,fitMethod,betaVar,modeVar)
        algorithm(obj);
        edgePointOrder = updateBICWeights(obj,maxDegreeBezier,maxCurvature,fitMethod,betaVar,modeVar);        
        stop = mergeClustersBIC(obj,maxDegreeBezier,maxCurvature,fitMethod,betaVar,modeVar);
        assignPointsToModels(obj,nSigmaThreshold);
        assignPointsToModels2(obj);
        dissolveClustersSmallerThan(obj,sizeThreshold);
        dissolveModelsLessDenseThan(obj,densityThreshold);
        dissolveModelsShorterThan(obj,thresholdLength);
        removeModels(obj,idx);
        saveClustersToHistory(obj);
        saveModelsToHistory(obj);
        saveEdgesToHistory(obj);
        idxPointsEdge = dissolveClustersCloseToEdge(obj,edgeDepth);
        filterGeomClusters(obj,filterData,nSigmaDist);
       
    end
    
    methods(Static = true)
        var = modelComponentVarianceWeighted(res,weights,varMin);
        var = modelDistanceVarianceWeightedWithPrior(res,weights,betaVar,modeVar);
        nParam = nParamModel(modelType,fitMethod);
        nParam = nParamModelMixture(modelType1,modelType2,fitMethod);
        logL = logLikelihoodModelAlongAway2D(res,weights,sigma,modelLength);
        logL = logLikelihoodModelAlongAway(res,weights,sigmaDistance,modelLength);
        logL = logLikelihoodModelMixture(logL1,logL2,n1,n2);
        bic = computeBIC(logL,nParam,nPoints);
    end
    
    methods (Access = private)
        weight = geomClusterWeight(obj,p1,p2,sigmaFilament,modelLength,angleThreshold);
    end
    
end

