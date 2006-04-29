function [residuals, Jacobian] = track_lsqnonlinFitFcn(parameters,sourceInfo,targetInfo,movieFrame,constants)
% track_lsqnonlinFitFcn is the fitting function for lsqnonlin

% the output is the catenated list of residuals, while the Jacobian is
% grad(img)*dW/dparms, where img is the image, W is the warping function
% and parms are the parameters (currently just translation warp)

nTags = length(sourceInfo);

% generate new targetInfo
[sourceInfo, targetInfo] = extractIntensities(sourceInfo, targetInfo,...
    movieFrame, constants, parameters, constants.gradientOption);

%---- debug
% if isfield(constants,'movieNoise')
% for iTag = 1:nTags
%     residuals = ...
%         targetInfo(iTag).deltaInt(targetInfo(iTag).goodIdx);
%     dof = length(residuals)-length(parameters);
%     % since we did a subtraction of two noisy signals, the
%     % noise should increase by sqrt(2)
%     % potential improvement: calculate robust second moment
%     sigmaResidual2(iTag) = (sum(residuals.^2)/(dof*2));
% 
%     % f-test
%     fRatio = sigmaResidual2(iTag)/constants.movieNoise(constants.currentTarget)^2;
%     fProb  = fcdf(fRatio,dof,numel(movieFrame)/2);
%     isSuccess(iTag) = fProb < 0.95;
% 
% end
% else
%     isSuccess = [NaN, NaN];
%     sigmaResidual2 = [NaN, NaN];
% end
% debug


% if 1% constants.verbose > 1
%     disp(sprintf(['parameters =',...
%         repmat(' %f',1,length(parameters(:))),...
%         ' - success %i %i sigmaResidual %f %f'],parameters',isSuccess,sigmaResidual2))
% end


% make list of residuals
switch constants.gaussPixOnly
    case 0
        residuals = cat(1,targetInfo.deltaInt);
    case 1
        % use pixels belonging to Gauss only
        residuals = cat(1,targetInfo.goodIdx);
        nResiduals = 0;
        for iTag = 1:nTags
            tagResiduals = targetInfo(iTag).deltaInt(sourceInfo(iTag).goodIdx);
            startIndex = nResiduals + 1;
            nResiduals = nResiduals + length(tagResiduals);
            residuals(startIndex:nResiduals) = tagResiduals;
        end
end

disp(sprintf(['parameters =',...
        repmat(' %f',1,length(parameters(:))),...
        '- res %f'],parameters',sum(residuals.^2)))



if any(isnan(residuals))
    'isnan residuals'
    keyboard
end


% Make the Jacobian to speed up convergence
if nargout > 1
    Jacobian = blkdiag(targetInfo.gradient);
    %Jacobian = []
end

%debug
if isfield(constants,'doSum') && constants.doSum==1
    residuals = sum(residuals.^2);
end

% more debug - needs gradient
% trackerDebugMakeMovie(sourceInfo,targetInfo);

