function updateModels(obj,maxCurvature,fitMethod,betaVar,modeVar)

%
% Updates the models of all clusters
%

% Copy data
obj_data_points = obj.data.points;
obj_data_error = obj.data.error;
obj_data_clusters = obj.data.clusters;
obj_data_modelType = obj.data.modelType;
obj_data_modelLength = obj.data.modelLength;
obj_data_modelRes = obj.data.modelRes;
obj_data_modelProj = obj.data.modelProj;
obj_data_modelBezCP = obj.data.modelBezCP;
obj_data_modelIsOutOfDate = obj.data.modelIsOutOfDate;
obj_data_modelVar = obj.data.modelVar;

parfor c=1:obj.data.nClusters
    if obj_data_modelIsOutOfDate(c)
        points = obj_data_points(obj_data_clusters{c},:);
        locPrec = obj_data_error(obj_data_clusters{c},:);
        
        % Reorder the points
        [~,~,~,idxProj,~] = projectPointsOntoLine(points,obj_data_modelBezCP{c}(1,:),obj_data_modelBezCP{c}(1,:)-obj_data_modelBezCP{c}(end,:));
        points = points(idxProj,:);
        locPrec = locPrec(idxProj,:);
        obj_data_clusters(c) = {obj_data_clusters{c}(idxProj)};
        
        % Update the model
        pointWeights = 1./locPrec;
        
        switch(fitMethod)
            case {1,2,5} % std. 3D, std. 2D
                [cP,t,res] = TLSFitBezierWeightedConstrainedCP(points,pointWeights,obj_data_modelType(c),maxCurvature);
                
            case {3,4} % snakes 3D, snakes 2D
                sigmaX = diag(locPrec(1,:).^2);
                [cP,t,res] = snakeBasedBezierFitUnderConstraint3(points,obj_data_modelType(c),modeVar*maxCurvature,'SigmaX',sigmaX);
                
            otherwise
                disp('Process: WARNING: Unknown fit method!');
        end
        
        % Compute the variance of the model
        var = Processor.modelDistanceVarianceWeightedWithPrior(res,pointWeights,betaVar,modeVar);
        
        [t,sortIdx] = sort(t);
        res = res(sortIdx,:);
        
        length = lengthBezier(cP);
        
        obj_data_modelLength(c) = length;
        obj_data_modelRes(c) = {res};
        obj_data_modelProj(c) = {t};
        obj_data_modelBezCP(c) = {cP};
        obj_data_modelVar(c) = var;
        obj_data_modelIsOutOfDate(c) = false;
        
        % Sort the points
        obj_data_clusters(c) = {obj_data_clusters{c}(sortIdx)};
    end % if
end % for

% Update data
obj.data.clusters = obj_data_clusters;
obj.data.modelLength = obj_data_modelLength;
obj.data.modelRes = obj_data_modelRes;
obj.data.modelProj = obj_data_modelProj;
obj.data.modelBezCP = obj_data_modelBezCP;
obj.data.modelVar = obj_data_modelVar;
obj.data.modelIsOutOfDate = obj_data_modelIsOutOfDate;

disp('Process: Models updated!');

end
