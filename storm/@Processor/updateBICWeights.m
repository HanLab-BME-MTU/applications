function edgePointOrder = updateBICWeights(obj,maxDegreeBezier,maxCurvature,fitMethod,betaVar,modeVar)

global stormTimer__; % :-D

% % % if maxCurvature < 0
% % %     maxCurvature = -maxCurvature;
% % %     penalizeSplitModel = true;
% % %     disp('Process: maxCurvature is negative!')
% % % else
% % %     penalizeSplitModel = false;
% % % end

% Remove edges containing a small cluster
idx = ismember(obj.data.edges(:,1),0);
obj.data.edges = obj.data.edges(~idx,:);
idx = ismember(obj.data.edges(:,2),0);
obj.data.edges = obj.data.edges(~idx,:);

nEdges = size(obj.data.edges,1);

% Initialize and copy data
obj_data_edges = obj.data.edges;
obj_data_points = obj.data.points;
obj_data_error = obj.data.error;
obj_data_clusters = obj.data.clusters;
obj_data_modelBezCP = obj.data.modelBezCP;

obj_data_modelLength = obj.data.modelLength;
obj_data_modelRes = obj.data.modelRes;
obj_data_modelType = obj.data.modelType;

obj_data_weights = zeros(nEdges,1);
obj_edgeModelType = zeros(nEdges,1);
obj_edgeModelLength = zeros(nEdges,1);
obj_edgeModelRes = cell(nEdges,1);
obj_edgeModelProj = cell(nEdges,1);
obj_edgeModelBezCP = cell(nEdges,1);
obj_edgeModelVar = zeros(nEdges,1);

obj_edgePointOrder = cell(nEdges,1);

clusterSize = cellfun(@numel,obj.data.clusters);
edgeClusterSize = sum(clusterSize(obj.data.edges),2);
effNPoint = sum(edgeClusterSize);
avNPoint = mean(edgeClusterSize);
data = [effNPoint,avNPoint,nEdges];

stormTimer__.start('BIC Parfor');
parfor e=1:nEdges
% for e=1:nEdges
    % Get the cluster index
    cluster1 = obj_data_edges(e,1);
    cluster2 = obj_data_edges(e,2);
    
    % Get the point data
    points1 = obj_data_points(obj_data_clusters{cluster1},:);
    points2 = obj_data_points(obj_data_clusters{cluster2},:);
    points = [points1;points2];
    
    nPoints1 = size(points1,1);
    nPoints2 = size(points2,1);
    nPoints = size(points,1);
    
    locPrec1 = obj_data_error(obj_data_clusters{cluster1},:);
    locPrec2 = obj_data_error(obj_data_clusters{cluster2},:);
    locPrec = [locPrec1;locPrec2];
    
    pointWeights1 = 1./locPrec1;
    pointWeights2 = 1./locPrec2;
    pointWeights = 1./locPrec;
    
    % Get the Bezier control points for both submodels
    cP1 = obj_data_modelBezCP{cluster1};
    cP2 = obj_data_modelBezCP{cluster2};
    
    % Find the two most distant submodels endpoints
    pEnd = [cP1(1,:);cP1(end,:);cP2(1,:);cP2(end,:)];
    distances = createDistanceMatrix(pEnd,pEnd);
    [rowVal rowIdx] = max(distances,[],1);
    [~,colIdx] = max(rowVal);
    linePoint = pEnd(colIdx,:);
    orientation = pEnd(colIdx,:)-pEnd(rowIdx(colIdx),:);
    
    % Project the points and sort them
    [~,~,~,idxProj,~] = projectPointsOntoLine(points,linePoint,orientation);
    points = points(idxProj,:);
    pointWeights = pointWeights(idxProj,:);
    
    type1 = obj_data_modelType(cluster1);
    type2 = obj_data_modelType(cluster2);
    
    modelDeg = (1:min(nPoints-2,maxDegreeBezier))';
    
    % Fit different models
    switch(fitMethod)
        case {1,2} % std. 3D or std. 2D
            fun = @(a) TLSFitBezierWeightedConstrainedCP(points,pointWeights,a,maxCurvature);
            [cPCurve,t,res] = arrayfun(fun,modelDeg,'UniformOutput',false);
            
        case {3,4} % snakes 3D or snakes 2D
            sigmaX = diag(locPrec(1,:).^2);
            fun = @(a) snakeBasedBezierFitUnderConstraint3(points,a,modeVar*maxCurvature,'SigmaX',sigmaX);
            [cPCurve,t,res] = arrayfun(fun,modelDeg,'UniformOutput',false);
            
        otherwise
            disp('BicCurveMatcher: WARNING: Unknown fit method!');
    end
    
    % Split model
    nParam12 = Processor.nParamModelMixture(type1,type2,fitMethod);
    
    length1 = obj_data_modelLength(cluster1);
    length2 = obj_data_modelLength(cluster2);
    
    res1 = obj_data_modelRes{cluster1};
    res2 = obj_data_modelRes{cluster2};

    var1 = Processor.modelDistanceVarianceWeightedWithPrior(res1,pointWeights1,betaVar,modeVar);
    var2 = Processor.modelDistanceVarianceWeightedWithPrior(res2,pointWeights2,betaVar,modeVar);
    
    logL1 = Processor.logLikelihoodModelAlongAway(res1,pointWeights1,sqrt(var1),length1);
    logL2 = Processor.logLikelihoodModelAlongAway(res2,pointWeights2,sqrt(var2),length2);
    logL12 = Processor.logLikelihoodModelMixture(logL1,logL2,nPoints1,nPoints2);
    
    % The merged models
    nParam = Processor.nParamModel(modelDeg,fitMethod);
    
    length = cellfun(@lengthBezier,cPCurve);
     
    var = cellfun(@(a) Processor.modelDistanceVarianceWeightedWithPrior(a,pointWeights,betaVar,modeVar),res); 

    logL = cellfun(@(a,b,c) Processor.logLikelihoodModelAlongAway(a,pointWeights,b,c),res,num2cell(sqrt(var)),num2cell(length));
    
    % % %     if penalizeSplitModel
    % % %         gapFactor = 1.25;
    % % %         length(maxDegreeBezier+1) = max(lengthC1+lengthC2,min((lengthC1+lengthC2)*gapFactor,max(length(1:maxDegreeBezier))));
    % % %     else
    % % %         length(maxDegreeBezier+1) = lengthC1+lengthC2;
    % % %     end
    % % %     % Penalize the split model
    % % %     logL(maxDegreeBezier+1) = logL(maxDegreeBezier+1) + nPoints*log(0.9);
    
    % Compute the BIC
    bic = Processor.computeBIC(logL,nParam,nPoints);
    bic12 = Processor.computeBIC(logL12,nParam12,nPoints);
        
    % Find the model with the smallest BIC
    [~,idxModel] = min([bic;bic12]);
                
    % If it is not the split model, then compute the
    % weight otherwise disable the edge

    if idxModel ~= modelDeg(end)+1       
        [t,idxFit] = sort(t{idxModel});
        obj_data_weights(e) = bic12 - bic(idxModel);
        obj_edgeModelType(e,1) = idxModel;
        obj_edgeModelLength(e,1) = length(idxModel);
        obj_edgeModelRes{e,1} = res{idxModel}(idxFit,:);
        obj_edgeModelProj{e,1} = t;
        obj_edgeModelBezCP{e,1} = cPCurve{idxModel};
        obj_edgeModelVar(e) = var(idxModel);
        obj_edgeModelVar(e) = 1;
        obj_edgePointOrder{e} = idxProj(idxFit);
    else
        obj_data_weights(e) = -1; % Disable the edge
    end
    
end % for
stormTimer__.stop('BIC Parfor',data);

obj.data.weights = obj_data_weights;
obj.edgeModelType = obj_edgeModelType;
obj.edgeModelLength = obj_edgeModelLength;
obj.edgeModelRes = obj_edgeModelRes;
obj.edgeModelProj = obj_edgeModelProj;
obj.edgeModelBezCP = obj_edgeModelBezCP;
obj.edgeModelVar = obj_edgeModelVar;
edgePointOrder = obj_edgePointOrder;

end % function
