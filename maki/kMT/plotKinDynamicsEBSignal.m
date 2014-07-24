function plotKinDynamicsEBSignal(sisterListEB,poleInfo,varargin)
%PLOTKINDYNAMICSEBSIGNAL makes plots of kinetochore dynamics and EB signal
%
% SYNOPSIS: plotKinDynamicsEBSignal(sisterListEB,poleInfo,varargin)
%
% INPUT sisterListEB: Output of detectkEBs.
%       poleInfo: Output of getSpindlePolesEB.
%       varargin: Optional input variables in the form of variable
%                 name/value pairs. Currently this includes:
%                   'timeLapse' - time between frames. Default: 1.
%                   'timeUnit' - time unit. Default: 'frames'.
%                   'pixelSize' - pixel size. Default: 1.
%                   'spaceUnit' - space unit. Default: 'pixels'.
%
% OUTPUT
%
% created by: Khuloud Jaqaman, June 2013

%% TEST INPUT & READ PARAMETERS

% Input check
assert( (isvector(sisterListEB) && isstruct(sisterListEB(1))) );
assert( (isvector(poleInfo) && isstruct(poleInfo(1))) );
ip = inputParser;
ip.addParamValue('timeLapse',1,@isscalar);
ip.addParamValue('timeUnit','frames',@ischar);
ip.addParamValue('pixelSize',1,@isscalar);
ip.addParamValue('spaceUnit','pixels',@ischar);
ip.parse(varargin{:});

% Get parameters
timeLapse = ip.Results.timeLapse;
timeUnit = ip.Results.timeUnit;
pixelSize = ip.Results.pixelSize;
spaceUnit = ip.Results.spaceUnit;

%number of sisters
numSister = length(sisterListEB);

%% BASIC KINETOCHORE DYNAMICS AND EB SIGNAL INFORMATION

%get number of frames and number of poles in spindle
numFrame = length(poleInfo);
numPole = size(poleInfo(1).xCoord,1);

%extract pole positions and distance between poles
coordPole = zeros(numFrame,3,numPole);
tmp = cat(2,poleInfo.xCoord);
coordPole(:,1,:) = tmp(:,1:2:end)';
tmp = cat(2,poleInfo.yCoord);
coordPole(:,2,:) = tmp(:,1:2:end)';
if isfield(poleInfo,'zCoord')
    tmp = cat(2,poleInfo.zCoord);
    coordPole(:,3,:) = tmp(:,1:2:end)';
end
% distPoles = sqrt( sum( (coordPole(:,:,1) - coordPole(:,:,2)).^2,2 ) );

%go over sisters and make plots
for iSis = 1 : numSister
    
    %get which pole each kinetochore is attached to
    sisterPoleAssign = sisterListEB(iSis).poleAssign12;
    
    %get sister coordinates
    coordsSis1 = sisterListEB(iSis).coords1;
    coordsSis2 = sisterListEB(iSis).coords2;
    
    %calculate sister-to-pole distances
    vecSis1PoleNear = coordsSis1(:,1:3) - coordPole(:,:,sisterPoleAssign(1,1));
    distSis1PoleNear = sqrt( sum( vecSis1PoleNear.^2,2 ) );
    sigmaDist = sqrt( sum( (vecSis1PoleNear.*coordsSis1(:,4:6)./repmat(distSis1PoleNear,1,3)).^2,2 ) );
    distSis1PoleNear = [distSis1PoleNear sigmaDist] * pixelSize;
    vecSis1PoleFar = coordsSis1(:,1:3) - coordPole(:,:,sisterPoleAssign(1,2));
    distSis1PoleFar = sqrt( sum( vecSis1PoleFar.^2,2 ) );
    sigmaDist = sqrt( sum( (vecSis1PoleFar.*coordsSis1(:,4:6)./repmat(distSis1PoleFar,1,3)).^2,2 ) );
    distSis1PoleFar = [distSis1PoleFar sigmaDist] * pixelSize;
    vecSis2PoleNear = coordsSis2(:,1:3) - coordPole(:,:,sisterPoleAssign(1,2));
    distSis2PoleNear = sqrt( sum( vecSis2PoleNear.^2,2 ) );
    sigmaDist = sqrt( sum( (vecSis2PoleNear.*coordsSis2(:,4:6)./repmat(distSis2PoleNear,1,3)).^2,2 ) );
    distSis2PoleNear = [distSis2PoleNear sigmaDist] * pixelSize;
    vecSis2PoleFar = coordsSis2(:,1:3) - coordPole(:,:,sisterPoleAssign(1,1));
    distSis2PoleFar = sqrt( sum( vecSis2PoleFar.^2,2 ) );
    sigmaDist = sqrt( sum( (vecSis2PoleFar.*coordsSis2(:,4:6)./repmat(distSis2PoleFar,1,3)).^2,2 ) );
    distSis2PoleFar = [distSis2PoleFar sigmaDist] * pixelSize;
    
    %get sister-to-sister distances
    distSis12 = sisterListEB(iSis).distances * pixelSize;
    
    %get EB signals
    signalEBSis1 = sisterListEB(iSis).kEBamp1(:,1);
    signalEBSis2 = sisterListEB(iSis).kEBamp2(:,1);
    
    %make plots
    h = figure('Name',['Sister Pair ' num2str(iSis)]);
    hold on
    
    %plot 1 - distances relative to pole that sister 1 is attached to
    subplot(2,3,1), hold on
    plot(distSis1PoleNear(:,1))
    myErrorbar(distSis1PoleNear(:,1),distSis1PoleNear(:,2))
    plot(distSis2PoleFar(:,1),'r')
    myErrorbar(distSis2PoleFar(:,1),distSis2PoleFar(:,2))
    legend('1','2')
    xlabel(['Time (' timeUnit ')'])
    ylabel(['Distance to Pole that Sister 1 is attached to (' spaceUnit ')'])
    
    %plot 2 - distances relative to pole that sister 2 is attached to
    subplot(2,3,4), hold on
    plot(distSis1PoleFar(:,1))
    myErrorbar(distSis1PoleFar(:,1),distSis1PoleFar(:,2))
    plot(distSis2PoleNear(:,1),'r')
    myErrorbar(distSis2PoleNear(:,1),distSis2PoleNear(:,2))
    legend('1','2')
    xlabel(['Time (' timeUnit ')'])
    ylabel(['Distance to Pole that Sister 2 is attached to (' spaceUnit ')'])
    
    %plot 3 - distances between each sister and its pole
    subplot(2,3,2), hold on
    plot(distSis1PoleNear(:,1))
    myErrorbar(distSis1PoleNear(:,1),distSis1PoleNear(:,2))
    plot(distSis2PoleNear(:,1),'r')
    myErrorbar(distSis2PoleNear(:,1),distSis2PoleNear(:,2))
    legend('1','2')
    xlabel(['Time (' timeUnit ')'])
    ylabel(['Distance between each sister and its pole (' spaceUnit ')'])
    
    %plot 4 - EB signal
    subplot(2,3,5), hold on
    plot(signalEBSis1)
    plot(signalEBSis2,'r')
    legend('1','2')
    xlabel(['Time (' timeUnit ')'])
    ylabel('EB signal intensity (a.u.)')
    
    %plot 5 - distance between sisters
    subplot(2,3,3), hold on
    plot(distSis12(:,1),'k')
    myErrorbar(distSis12(:,1),distSis12(:,2))
    xlabel(['Time (' timeUnit ')'])
    ylabel(['Distance between sisters (' spaceUnit ')'])
    
end

