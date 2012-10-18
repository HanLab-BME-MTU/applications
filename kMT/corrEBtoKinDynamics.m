function measurementsEB = corrEBtoKinDynamics(sisterListEB,poleInfo,varargin)
%CORREBTOKINDYNAMICS correlates EB signal at kinetochores to kinetochore dynamics
%
% SYNOPSIS: corrEBtoKinDynamics(sisterListEB,poleInfo,spindleAxis,varargin)
%
% INPUT sisterListEB : Output of detectkEBs.
%       poleInfo: Structure array in the form of movieInfo (output of
%                 detectSubResFeatures2D_StandAlone) with spindle pole
%                 position and amplitude information. Currently, output of
%                 getSpindlePolesEB.
%       varargin     : Optional input variables in the form of variable
%                      name/value pairs. Currently this includes:
%                           ...
%
% OUTPUT 
%
% created by: Khuloud Jaqaman, October 2012

%% TEST INPUT & READ PARAMETERS

% Input check
assert(isvector(sisterListEB) && isstruct(sisterListEB(1)));
assert(isvector(poleInfo) && isstruct(poleInfo(1)));
ip = inputParser;
% ip.addParamValue('radiusEB',7,@isscalar);
ip.parse(varargin{:});

% % Get parameters
% radiusEB = ip.Results.radiusEB;

%number of sisters
numSister = length(sisterListEB);

%% BASIC KINETOCHORE DYNAMICS AND EB SIGNAL INFORMATION

%get number of frames and number of poles in spindle
numFrame = length(poleInfo);
numPole = size(poleInfo(1).xCoord,1);

%extract pole positions
coordPole = zeros(numFrame,3,numPole);
tmp = cat(2,poleInfo.xCoord);
coordPole(:,1,:) = tmp(:,1:2:end)';
tmp = cat(2,poleInfo.yCoord);
coordPole(:,2,:) = tmp(:,1:2:end)';
if isfield(poleInfo,'zCoord')
    tmp = cat(2,poleInfo.zCoord);
    coordPole(:,3,:) = tmp(:,1:2:end)'; %#ok<NASGU>
end

%retrieve which pole each kinetochore is attached to
sisterPoleAssign = NaN(numFrame,2,numSister);
for iSis = 1 : numSister
    sisterPoleAssign(:,:,iSis) = sisterListEB(iSis).poleAssign12;
end

%various measurements characterising kinetochore dynamics
dispSisPair = NaN(numFrame-1,3,2,numSister); %frame-to-frame displacement
dispMagSisPair = NaN(numFrame-1,2,numSister); %displacement magnitude
vecSisPair = NaN(numFrame,3,2,numSister); %vector connecting kinetochore to its pole
distSisPair = NaN(numFrame,2,numSister); %magnitude of above vector
dispProjSisPair = NaN(numFrame-1,2,numSister); %projection of frame-to-frame displacement onto vector connecting kinetochore to pole
for iSis = 1 : numSister
    for j = 1 : 2
        eval(['coordsSisJ = sisterListEB(iSis).coords' num2str(j) '(:,1:3);'])
        dispSisPair(:,:,j,iSis) = diff(coordsSisJ,[],1);
        dispMagSisPair(:,j,iSis) = normList(dispSisPair(:,:,j,iSis));
        eval(['vecSisPair(:,:,j,iSis) = sisterListEB(iSis).vecFromPole' num2str(j) '(:,1:3);'])
        distSisPair(:,j,iSis) = normList(vecSisPair(:,:,j,iSis));
        dispProjSisPair(:,j,iSis) = sum(dispSisPair(:,:,j,iSis) .* vecSisPair(1:end-1,:,j,iSis),2) ./ distSisPair(1:end-1,j,iSis);
    end
end

%put together EB signals to facilitate further analysis
ebSignalPair = NaN(numFrame,2,numSister);
for iSis = 1 : numSister
    for j = 1 : 2
        eval(['ebSignalPair(:,j,iSis) = sisterListEB(iSis).kEBamp' num2str(j) '(:,1);'])
    end
end

%% CORRELATE EB SIGNAL TO KINETOCHORE DYNAMICS

minDisp = 0.5;
minDispAway = minDisp;
minDispToward = -minDisp;

%get EB signals away and toward
ebSignalForDisp = ebSignalPair(1:end-1,:,:);
[ebSignalAway,ebSignalToward] = deal(NaN(numFrame-1,2,numSister));
ebSignalAway(dispProjSisPair>minDispAway) = ebSignalForDisp(dispProjSisPair>minDispAway);
ebSignalToward(dispProjSisPair<minDispToward) = ebSignalForDisp(dispProjSisPair<minDispToward);

%probability of having an EB comet if moving away from or toward pole, per
%kinetochore

numCometAway = squeeze(sum(ebSignalAway>0,1));
numCometToward = squeeze(sum(ebSignalToward>0,1));
numTimesMoveAway = squeeze(sum(dispProjSisPair>minDispAway,1));
numTimesMoveToward = squeeze(sum(dispProjSisPair<minDispToward,1));

probEBCometAway = numCometAway ./ numTimesMoveAway;
probEBCometToward = numCometToward ./ numTimesMoveToward;

%EB ratio within each kinetochore between moving away from and moving
%toward pole

meanEBSignalAway = squeeze(nanmean(ebSignalAway,1));
meanEBSignalToward = squeeze(nanmean(ebSignalToward,1));

ratioEBAwayToward = meanEBSignalAway ./ meanEBSignalToward;

% %cross-correlation between kinetochore displacement and EB signal
% 
% %collect time series in the appropriate format
% dispStruct = repmat(struct('observations',[]),2,numSister);
% ebStruct = repmat(struct('observations',[]),2,numSister);
% for iSis = 1 : numSister
%     for j = 1 : 2
%         dispStruct(j,iSis).observations = dispProjSisPair(:,j,iSis);
%         ebStruct(j,iSis).observations = ebSignalForDisp(:,j,iSis);
%     end
% end
% dispStruct = dispStruct(:);
% ebStruct = ebStruct(:);
% 
% crossCorrDispEB = crossCorr(dispStruct,ebStruct,20,1);

%% OUTPUT

measurementsEB.numTimesKinMoveAway = numTimesMoveAway;
measurementsEB.numTimesKinMoveToward = numTimesMoveToward;
measurementsEB.probEBCometAway = probEBCometAway;
measurementsEB.probEBCometToward = probEBCometToward;
measurementsEB.ratioEBSignalAwayToward = ratioEBAwayToward;


