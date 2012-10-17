function [ratioEBAwayToward,crossCorrDispEB] = corrEBtoKinDynamics(sisterListEB,poleInfo,varargin)
%CORREBTOKINDYNAMICS correlates EB signal at kinetochores to kinetochore dynamics
%
% SYNOPSIS: corrEBtoKinDynamics(sisterListEB,poleInfo,spindleAxis,varargin)
%
% INPUT sisterListEB : Output of detectkEBs.
%       poleInfo: Structure array in the form of movieInfo (output of
%                 detectSubResFeatures2D_StandAlone) with spindle pole
%                 position and amplitude information. Currently, output of
%                 getSpindleAxisEB, but that is not essential.
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
    coordPole(:,3,:) = tmp(:,1:2:end)';
end 

%determine which spindle pole each kinetochore in a sister pair is attached to
%this assignment will depend on the spindle geometry (number of poles)
%currently code handles only bipolar spindles and a simplified monopolar
%case where all attachments are monotelic
%this assignment is simplistic, it could make use of the EB signal to make
%a more sophisticated assignment, even in bipolar spindles
%but this definitely works in cases of sisters aligned at the metaphase
%plate
sisterPoleAssign = NaN(numSister,2); %NaN will remain to mean not attached to any spindle pole
switch numPole
    
    case 1 %monopolar spindle
        
        for iSis = 1 : numSister
            distSis1 = nanmean(normList(sisterListEB(iSis).coords1(:,1:3)-coordPole(:,:,1)));
            distSis2 = nanmean(normList(sisterListEB(iSis).coords2(:,1:3)-coordPole(:,:,1)));
            if distSis2 < distSis1
                sisterPoleAssign(iSis,2) = 1;
            else
                sisterPoleAssign(iSis,1) = 1;
            end
                
        end
        
    case 2 %bipolar spindle
        
        for iSis = 1 : numSister
            distSis1 = nanmean(normList(sisterListEB(iSis).coords1(:,1:3)-coordPole(:,:,1)));
            distSis2 = nanmean(normList(sisterListEB(iSis).coords2(:,1:3)-coordPole(:,:,1)));
            if distSis2 < distSis1
                sisterPoleAssign(iSis,:) = [2 1];
            else
                sisterPoleAssign(iSis,:) = [1 2];
            end
        end
        
    otherwise %multipolar spindle
        
end

%various measurements characterising kinetochore dynamics
dispSisPair = NaN(numFrame-1,3,2,numSister); %frame-to-frame displacement
dispMagSisPair = NaN(numFrame-1,2,numSister); %displacement magnitude
vecSisPair = NaN(numFrame,3,2,numSister); %vector connecting kinetochore to its pole
distSisPair = NaN(numFrame,2,numSister); %magnitude of above vector
dispProjSisPair = NaN(numFrame-1,2,numSister); %projection of frame-to-frame displacement onto vector connecting kinetochore to pole
for iSis = 1 : numSister
    for j = 1 : 2
        poleSisJ = sisterPoleAssign(iSis,j);
        eval(['coordsSisJ = sisterListEB(iSis).coords' num2str(j) '(:,1:3);'])
        dispSisPair(:,:,j,iSis) = diff(coordsSisJ,[],1);
        dispMagSisPair(:,j,iSis) = normList(dispSisPair(:,:,j,iSis));
        if ~isnan(poleSisJ)
            vecSisPair(:,:,j,iSis) = coordsSisJ - coordPole(:,:,poleSisJ);
        end
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

%average EB signal ratio between moving away from and moving toward the pole

%get EB signals away and toward
[ebSignalAway,ebSignalToward] = deal(NaN(numFrame-1,2,numSister));
ebSignalTmp = ebSignalPair(1:end-1,:,:);
ebSignalAway(dispProjSisPair>0) = ebSignalTmp(dispProjSisPair>0);
ebSignalToward(dispProjSisPair<0) = ebSignalTmp(dispProjSisPair<0);

%calculate average signal per kinetochore
meanEBSignalAway = squeeze(nanmean(ebSignalAway,1));
meanEBSignalToward = squeeze(nanmean(ebSignalToward,1));

%calculate ratio per kinetochore
ratioEBAwayToward = meanEBSignalAway ./ meanEBSignalToward;

%cross-correlation between kinetochore displacement and EB signal

%collect time series in the appropriate format
dispStruct = repmat(struct('observations',[]),2,numSister);
ebStruct = repmat(struct('observations',[]),2,numSister);
for iSis = 1 : numSister
    for j = 1 : 2
        dispStruct(j,iSis).observations = dispProjSisPair(:,j,iSis);
        ebStruct(j,iSis).observations = ebSignalTmp(:,j,iSis);
    end
end
dispStruct = dispStruct(:);
ebStruct = ebStruct(:);

crossCorrDispEB = crossCorr(dispStruct,ebStruct,20,1);

%% OUTPUT

