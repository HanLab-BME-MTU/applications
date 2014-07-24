function initModels(obj,betaVar,modeVar)

obj_data_modelType = zeros(obj.data.nClusters,1);
obj_data_modelLength = zeros(obj.data.nClusters,1);
obj_data_modelRes = cell(obj.data.nClusters,1);
obj_data_modelProj = cell(obj.data.nClusters,1);
obj_data_modelBezCP = cell(obj.data.nClusters,1);
obj_data_modelVar = zeros(obj.data.nClusters,1);
obj.data.modelIsOutOfDate = false(obj.data.nClusters,1);

obj_data_points = obj.data.points;
obj_data_clusters = obj.data.clusters;
obj_data_error = obj.data.error;

parfor c=1:obj.data.nClusters
    % Compute the linear model without weights to determine the point order
    points = obj_data_points(obj_data_clusters{c},:);
    
    meanW = mean(1./obj_data_error(obj_data_clusters{c},:),1);
    [x0,a,~,~] = ls3dline(bsxfun(@times,points,meanW));
    
    % Determine the point order
    [~,~,~,idx,~] = projectPointsOntoLine(points,x0'./meanW,a'./meanW);
    
    % Sort the points
    obj_data_clusters(c) = {obj_data_clusters{c}(idx)};
    
    % Compute the linear model with weights
    points = obj_data_points(obj_data_clusters{c},:);
    pointWeights = 1./obj_data_error(obj_data_clusters{c},:);
    maxCurvature = 100; % Dummy value: Linear model
    [cPCurve,t,res] = TLSFitBezierWeightedConstrainedCP(points,pointWeights,1,maxCurvature);
    
    % Compute the variance of the model
    var = Processor.modelDistanceVarianceWeightedWithPrior(res,pointWeights,betaVar,modeVar);
    
    % Sort the points
    [t,sortIdx] = sort(t);
    res = res(sortIdx,:);
    
    % Compute the model length
    length = lengthBezier(cPCurve);
    
    obj_data_modelType(c) = 1;
    obj_data_modelLength(c) = length;
    obj_data_modelRes(c) = {res};
    obj_data_modelProj(c) = {t};
    obj_data_modelBezCP(c) = {cPCurve};
    obj_data_modelVar(c) = var;
    
    % Sort the points
    obj_data_clusters(c) = {obj_data_clusters{c}(sortIdx)};
end

obj.data.modelType = obj_data_modelType;
obj.data.modelLength = obj_data_modelLength;
obj.data.modelRes = obj_data_modelRes;
obj.data.modelProj = obj_data_modelProj;
obj.data.modelBezCP = obj_data_modelBezCP;
obj.data.modelVar = obj_data_modelVar;
obj.data.clusters = obj_data_clusters;

disp('Process: Models initialized!');

end
