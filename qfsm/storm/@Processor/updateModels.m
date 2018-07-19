function updateModels(obj,maxCurvature,fitMethod,betaVar,modeVar)
% function updateModels(obj,maxCurvature,fitMethod,betaVar,modeVar)
% SYNOPSIS:
% Updates the Bézier curve models of all clusters. 
%
% REQUIRED INPUTS:         
% - maxCurvature
% Depending on fitMethod this is the beta parameter of the snakes based fit
% or the maximum curvature for the standard fit.
%
% - fitMethod
% Different curve fitting methods can be specified.
% Standard 3D: 1
% Standard 2D: 2
% Snakes 3D: 3
% Snakes 2D: 4
%
% - betaVar
% The shape parameter of the variance prior.
%
% - modeVar
% The mode parameter of the variance prior.
% 
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES: 
% - obj.data.points
% - obj.data.clusters
% - obj.data.error
% - obj.data.modelBezCP
% - obj.data.modelType
% - obj.data.modelLength
% - obj.data.modelRes
% - obj.data.modelProj
% - obj.data.modelIsOutOfDate
% - obj.data.modelVar
%
% MODIFIED PROPERTIES: 
% - obj.data.clusters
% - obj.data.modelLength
% - obj.data.modelRes
% - obj.data.modelProj
% - obj.data.modelBezCP
% - obj.data.modelVar
% - obj.data.modelIsOutOfDate
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

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
        
        % Order the points
        [~,~,~,idxProj,~] = projectPointsOntoLine(points,obj_data_modelBezCP{c}(1,:),obj_data_modelBezCP{c}(1,:)-obj_data_modelBezCP{c}(end,:));
        points = points(idxProj,:);
        locPrec = locPrec(idxProj,:);
        obj_data_clusters(c) = {obj_data_clusters{c}(idxProj)};
        
        % Update the model
        pointWeights = 1./locPrec;
        
        switch(fitMethod)
            case {1,2} % std. 3D, std. 2D
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
