function sumOfSquaredResiduals = track_tomlabFitFcn(parameters,problem)
% track_tomlabFitFcn is the fitting function for tomlab

% read problem parameters
sourceInfo = problem.user.sourceInfo;
targetInfo = problem.user.targetInfo;
movieFrame = problem.user.movieFrame;
constants = problem.user.constants;

% generate new targetInfo
[sourceInfo, targetInfo] = extractIntensities(sourceInfo, targetInfo,...
    movieFrame, constants, parameters, constants.gradientOption);
if constants.verbose > 1
    disp(sprintf(['parameters =',...
        repmat(' %f',1,length(parameters(:)))],parameters'))
end

% make list of residuals
switch constants.gaussPixOnly
    case 0
        residuals = cat(1,targetInfo.deltaInt);
    case 1
        % use pixels belonging to Gauss only
        residuals = cat(1,targetInfo.goodIdx);
        nResiduals = 0;
        for iTag = 1:length(sourceInfo)
            tagResiduals = targetInfo(iTag).deltaInt(sourceInfo(iTag).goodIdx);
            startIndex = nResiduals + 1;
            nResiduals = nResiduals + length(tagResiduals);
            residuals(startIndex:nResiduals) = tagResiduals;
        end
end


sumOfSquaredResiduals = sum(residuals.^2);

