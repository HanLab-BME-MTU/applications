function [trackLft,trajClass,diffCoefGen,confRadAll,trajDiffMode,trajDiffCoef2,numDiffMode,...
    angleWithProtTmp,f2fDispTmp,paraDirDispTmp,perpDirDispTmp,paraProtDispTmp,...
    perpProtDispTmp,asymParamTmp,ratioDispDirTmp,ratioDispProtTmp] = ...
    collectParticleBehavior(tracksFinal,diffAnalysisRes,diffModeAnRes,directTrackChar)

%get the lifetime of each track segment
trackLft = getTrackSEL(tracksFinal,1);
trackLft = trackLft(:,3);

%get the number of track segments
numSegments = length(trackLft);

%% FROM ASYMMETRY AND DIFFUSION ANALYSIS ...

%get trajectory classifications
trajClassDiff = vertcat(diffAnalysisRes.classification);
trajClassAsym = trajClassDiff(:,1); %asymmetry classification
trajClassDiff = trajClassDiff(:,2); %diffusion classification

%WARNING
%MAKE trajClassAsym ALL NaN, SO THAT THERE IS NO ASYMMETRY CLASSIFICATION
trajClassAsym(:) = NaN;

%process classifications to make categories
if all(isnan(trajClassAsym))
    trajClass = trajClassDiff;
else
    trajClass = NaN(numSegments,1);
    trajClass(trajClassAsym==0&trajClassDiff==1) = 1; %isotropic+confined
    trajClass(trajClassAsym==0&trajClassDiff==2) = 2; %isotropic+free
    trajClass(trajClassAsym==0&trajClassDiff==3) = 3; %isotropic+directed
    trajClass(trajClassAsym==0&isnan(trajClassDiff)) = 4; %isotropic+undetermined
    trajClass(trajClassAsym==1) = 5; %linear+anything
end

%get diffusion coefficients
diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');

%get confinement radii
confRadAll = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');

%% FROM DIFFUSION MODE ANALYSIS ...

%get trajectory diffusion modes
trajDiffMode = vertcat(diffModeAnRes.diffMode);

%get trajectory diffusion coefficient
trajDiffCoef2 = vertcat(diffModeAnRes.diffCoef);

%get number of diffusion modes
numDiffMode = max(trajDiffMode);

%% FROM TRACKS DIRECTLY ...

%get direction of motion, angle with protrusion vector and various
%frame-to-frame displacement measures
angleWithProtTmp = vertcat(directTrackChar.angleWithProt);
f2fDispTmp       = vertcat(directTrackChar.f2fDisp);
paraDirDispTmp   = vertcat(directTrackChar.paraDirDisp);
perpDirDispTmp   = vertcat(directTrackChar.perpDirDisp);
paraProtDispTmp  = vertcat(directTrackChar.paraProtDisp);
perpProtDispTmp  = vertcat(directTrackChar.perpProtDisp);
asymParamTmp     = vertcat(directTrackChar.asymParam);

%calculate ratio of parallel to perpendicular displacements
ratioDispDirTmp = abs( paraDirDispTmp ./ perpDirDispTmp );
ratioDispProtTmp = abs( paraProtDispTmp ./ perpProtDispTmp );

