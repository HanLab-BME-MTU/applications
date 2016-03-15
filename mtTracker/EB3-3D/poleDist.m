function [dist,sphCoord,poleId,inliers,originProb,minProb,bestSphCoord,movieInfoSpindle] = poleDist(poleMovieInfo,pMovieInfo,varargin)
% The angle is defined centered on the closest pole. The base is either the
% 'absolute' referential of build around the axis between the closest pole and the second closest pole.
% WARNING:  designed and implemented for 2 poles, extensible to 3 poles.
% WARNING:  3D only
% 
% Philippe Roudot 2015
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('poleMovieInfo');
ip.addRequired('pMovieInfo');
ip.addParamValue('processFrames',[], @isnumeric);
ip.addParamValue('angleRef','poles', @ischar);
ip.addParamValue('anisotropy',[1 1 1], @isnumeric);
ip.addParamValue('Alpha',0.05, @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.parse(poleMovieInfo,pMovieInfo, varargin{:});

processFrames=[];
if isempty(ip.Results.processFrames)
    processFrames=1:length(poleMovieInfo);
else
    processFrames=ip.Results.processFrames;
end

anis=ip.Results.anisotropy;

dist=cell(1,numel(processFrames));
originProb=cell(1,numel(processFrames));
minProb=cell(1,numel(processFrames));
azimuth=cell(1,numel(processFrames));
elevation=cell(1,numel(processFrames));
rho=cell(1,numel(processFrames));
bazimuth=cell(1,numel(processFrames));
belevation=cell(1,numel(processFrames));
brho=cell(1,numel(processFrames));
poleId=cell(1,numel(processFrames));
inliers=cell(1,numel(processFrames));

movieInfoSpindle=pMovieInfo;
% Reference axis use for registration.
refInterPolarAxis=[ anis(1)*poleMovieInfo(1).xCoord(2,1)-anis(1)*poleMovieInfo(1).xCoord(1,1) ...
                    anis(2)*poleMovieInfo(1).yCoord(2,1)-anis(2)*poleMovieInfo(1).yCoord(1,1) ...
                    anis(3)*poleMovieInfo(1).zCoord(2,1)-anis(3)*poleMovieInfo(1).zCoord(1,1)];

refInterPolarAxisOrigin=[anis(1)*poleMovieInfo(1).xCoord(1,1) ...
                         anis(2)*poleMovieInfo(1).yCoord(1,1) ...
                         anis(3)*poleMovieInfo(1).zCoord(1,1)];                
                
% associated basis.
refVZ=refInterPolarAxis./repmat(sum(refInterPolarAxis.^2,2).^0.5,1,3);
refVX=[0*refVZ(1),refVZ(3),-refVZ(2)];
refVX=refVX./(sum(refVX.^2,2).^0.5);
refVY=cross(refVX,refVZ);


parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    
    poleNumber=size(poleMovieInfo(frameIdx).xCoord,1);
    particleNumber=size(pMovieInfo(frameIdx).xCoord,1);
    
    % Recenter each cart. coord. wrt each pole
    recenterCoord=zeros(particleNumber,poleNumber,3);
    recenterCoord(:,:,1)=anis(1)*(repmat(pMovieInfo(frameIdx).xCoord(:,1),1,poleNumber) - ...
                         repmat(poleMovieInfo(frameIdx).xCoord(:,1)',particleNumber,1));
    recenterCoord(:,:,2)=anis(2)*(repmat(pMovieInfo(frameIdx).yCoord(:,1),1,poleNumber) - ...
                         repmat(poleMovieInfo(frameIdx).yCoord(:,1)',particleNumber,1));
    recenterCoord(:,:,3)=anis(3)*(repmat(pMovieInfo(frameIdx).zCoord(:,1),1,poleNumber) - ...
                         repmat(poleMovieInfo(frameIdx).zCoord(:,1)',particleNumber,1));
    %% Compute distances
    dist{frameIdx}= sum(recenterCoord.^2,3).^0.5;  
    
    %% compute the closest pole ( designed for N poles)
    originProb{frameIdx}=dist{frameIdx}./repmat(sum(dist{frameIdx},2),1,poleNumber);
    minProb{frameIdx}=min(originProb{frameIdx},[],2);   
    [sortMinDist,sortedPoleIdx]=sort(dist{frameIdx},2);
    minDist=sortMinDist(:,1);
    closestPoleIdx=sortedPoleIdx(:,1);
    poleId{frameIdx}=closestPoleIdx; 
    
    Idx=sub2ind(size(recenterCoord),[1:particleNumber 1:particleNumber 1:particleNumber],[closestPoleIdx' closestPoleIdx' closestPoleIdx'],[ones(1,particleNumber) 2*ones(1,particleNumber) 3*ones(1,particleNumber)]);
    recenterCloserCoord=zeros(particleNumber,3);
    recenterCloserCoord(:)=recenterCoord(Idx);
                     
    % compute the azimuth, elevation and rho for each poles. 
    %% Two bases, one for each poles. 
    P1P2=zeros(1,3);
    P1P2(1)=anis(1)*(poleMovieInfo(frameIdx).xCoord(2,1)-poleMovieInfo(frameIdx).xCoord(1,1));
    P1P2(2)=anis(2)*(poleMovieInfo(frameIdx).yCoord(2,1)-poleMovieInfo(frameIdx).yCoord(1,1));
    P1P2(3)=anis(3)*(poleMovieInfo(frameIdx).zCoord(2,1)-poleMovieInfo(frameIdx).zCoord(1,1));

    vZ1=P1P2./repmat(sum(P1P2.^2,2).^0.5,1,3);
    vX1=[0*vZ1(:,1),vZ1(:,3),-vZ1(:,2)];
    vX1=vX1./repmat(sum(vX1.^2,2).^0.5,1,3);
    vY1=cross(vX1,vZ1);
    
    vZ2=-vZ1;
    vX2=[0*vZ2(:,1),vZ2(:,3),-vZ2(:,2)];
    vX2=vX2./repmat(sum(vX2.^2,2).^0.5,1,3);
    vY2=cross(vX2,vZ2);
    
    %% transform
    transformCoordPole=zeros(size(recenterCoord));
    for i=1:particleNumber
       transformCoordPole(i,1,:)=squeeze(recenterCoord(i,1,:))'*[vX1' vY1' vZ1'];
       transformCoordPole(i,2,:)=squeeze(recenterCoord(i,2,:))'*[vX2' vY2' vZ2'];
    end
    
    transformCoordPoleBestPoles=zeros(particleNumber,3);
    transformCoordPoleBestPoles(:)=transformCoordPole(Idx);
    
    [azimuth1,elevation1,rho1]=cart2sph(transformCoordPole(:,1,1),transformCoordPole(:,1,2),transformCoordPole(:,1,3));
    [azimuth2,elevation2,rho2]=cart2sph(transformCoordPole(:,2,1),transformCoordPole(:,2,2),transformCoordPole(:,2,3));
    [azimuthBest,elevationBest,rhoBest]=cart2sph(transformCoordPoleBestPoles(:,1),transformCoordPoleBestPoles(:,2),transformCoordPoleBestPoles(:,3));
    
    azimuth{frameIdx}=[azimuth1 azimuth2];
    elevation{frameIdx}=[elevation1 elevation2];
    rho{frameIdx}=[rho1 rho2];
    
    bazimuth{frameIdx}=azimuthBest;
    belevation{frameIdx}=elevationBest;
    brho{frameIdx}=rhoBest;
    
    %% Outlier estimation, flag particle that are too far from theire 
    %  associated pole through an otsu estimation to separate both
    %  population( works best when the outlier population is large).
    % If few outlier, switch to population analysis (see below), might need
    % an adaptive switch here. 
    %     firstPoleMinDist=minDist(closestPoleIdx==1);
    %     medDev=firstPoleMinDist-median(firstPoleMinDist);
    %     globalThreshold=median(firstPoleMinDist)+3*1.4628*median(medDev(medDev>0));
    
    l=graythresh(double(minDist)/max(minDist))*max(minDist);
    inliers{frameIdx}=minDist<l;
      
    %% Registration with the first frame spindle orientation
    % Reference use for registration.
    frameRefInterPolarAxis=[anis(1)*poleMovieInfo(frameIdx).xCoord(2,1)-anis(1)*poleMovieInfo(frameIdx).xCoord(1,1) ...
        anis(2)*poleMovieInfo(frameIdx).yCoord(2,1)-anis(2)*poleMovieInfo(frameIdx).yCoord(1,1) ...
        anis(3)*poleMovieInfo(frameIdx).zCoord(2,1)-anis(3)*poleMovieInfo(frameIdx).zCoord(1,1)];
    
    frameRefVZ=frameRefInterPolarAxis./repmat(sum(frameRefInterPolarAxis.^2,2).^0.5,1,3);
    frameRefVX=[0*frameRefVZ(1),frameRefVZ(3),-frameRefVZ(2)];
    frameRefVX=frameRefVX./repmat(sum(frameRefVX.^2,2).^0.5,1,3);
    frameRefVY=cross(frameRefVX,frameRefVZ);
    
    refBaseSwitchMatrix=[frameRefVX' frameRefVY' frameRefVZ']*[refVX' refVY' refVZ']^(-1);
    
    % Project on both adaptive basis
    transformCoordAxisRef=zeros(size(recenterCloserCoord));
    for i=1:particleNumber
        transformCoordAxisRef(i,:)=refInterPolarAxisOrigin+squeeze(recenterCoord(i,1,:))'*refBaseSwitchMatrix;
    end
    
    movieInfoSpindle(frameIdx).xCoord=[transformCoordAxisRef(:,1)/anis(1)  pMovieInfo(frameIdx).xCoord(:,2)];
    movieInfoSpindle(frameIdx).yCoord=[transformCoordAxisRef(:,2)/anis(2)  pMovieInfo(frameIdx).yCoord(:,2)];
    movieInfoSpindle(frameIdx).zCoord=[transformCoordAxisRef(:,3)/anis(3)  pMovieInfo(frameIdx).zCoord(:,2)];   
end

sphCoord.azimuth=azimuth;
sphCoord.elevation=elevation;
sphCoord.rho=rho;

bestSphCoord.azimuth=bazimuth;
bestSphCoord.elevation=belevation;
bestSphCoord.rho=brho;





