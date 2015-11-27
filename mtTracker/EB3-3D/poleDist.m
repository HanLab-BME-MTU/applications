function [dist,poleId,inliers,originProb,minProb,azimuth,elevation,rho,movieInfoSpindle] = poleDist(poleMovieInfo,pMovieInfo,varargin)
% The angle is defined centered on the closest pole. The base is either the
% 'absolute' referential of build around the axis between the closest pole and the second closest pole.
% WARNING:  designed and implemented for N poles, but tested with 2 poled
%           only. 
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
    
    % Compute distances
    dist{frameIdx}= sum(recenterCoord.^2,3).^0.5;  
    
    % compute the closest pole
    originProb{frameIdx}=dist{frameIdx}./repmat(sum(dist{frameIdx},2),1,poleNumber);
    minProb{frameIdx}=min(originProb{frameIdx},[],2);   
    [minDist,closestPoleIdx]=min(dist{frameIdx},[],2);
    [sortMinDist,sortedPoleIdx]=sort(dist{frameIdx},2);
    minDist=sortMinDist(:,1);
    closestPoleIdx=sortedPoleIdx(:,1);
    poleId{frameIdx}=closestPoleIdx;
    
    
%    reshapCloserPoleIdx=((repmat(minDist,[1,2,3])==repmat(dist{frameIdx},[1,1,3])))
%    recenterCloserCoord=zeros(particleNumber,3);
%    recenterCloserCoord(:)=recenterCoord(reshapCloserPoleIdx);
    Idx=sub2ind(size(recenterCoord),[1:particleNumber 1:particleNumber 1:particleNumber],[closestPoleIdx' closestPoleIdx' closestPoleIdx'],[ones(1,particleNumber) 2*ones(1,particleNumber) 3*ones(1,particleNumber)]);
    recenterCloserCoord=zeros(particleNumber,3);
    recenterCloserCoord(:)=recenterCoord(Idx);
    
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
    
    if(strcmp(ip.Results.angleRef,'poles'))
        %% Define the an interpolar axis associated to each particle
        poleCoord=zeros(particleNumber,poleNumber,3);
        poleCoord(:,:,1)=repmat(anis(1)*(poleMovieInfo(frameIdx).xCoord(:,1)'),particleNumber,1);
        poleCoord(:,:,2)=repmat(anis(2)*(poleMovieInfo(frameIdx).yCoord(:,1)'),particleNumber,1);
        poleCoord(:,:,3)=repmat(anis(3)*(poleMovieInfo(frameIdx).zCoord(:,1)'),particleNumber,1);
        
        bestPoles=zeros(particleNumber,3);
        bestPoles(:)=poleCoord(Idx);
        
        secClosestPoleIdx=sortedPoleIdx(:,2);
        Idx=sub2ind(size(recenterCoord),[1:particleNumber 1:particleNumber 1:particleNumber],[secClosestPoleIdx' secClosestPoleIdx' secClosestPoleIdx'],[ones(1,particleNumber) 2*ones(1,particleNumber) 3*ones(1,particleNumber)]);
        
        secondBestPoles=zeros(particleNumber,3);
        secondBestPoles(:)=poleCoord(Idx);
        
        interPolarAxis=(secondBestPoles-bestPoles);
        
        %% change coordinate for each EB3 (Z is along the polar axis)
        % Create a New Base
        vZ=interPolarAxis./repmat(sum(interPolarAxis.^2,2).^0.5,1,3);
        vX=[0*vZ(:,1),vZ(:,3),-vZ(:,2)];
        vX=vX./repmat(sum(vX.^2,2).^0.5,1,3);
        vY=cross(vX,vZ);
      
        
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
        transformCoordPoleRef=zeros(size(recenterCloserCoord));
        transformCoordAxisRef=zeros(size(recenterCloserCoord));
        for i=1:particleNumber
            transformCoordPoleRef(i,:)=recenterCloserCoord(i,:)*[vX(i,:)' vY(i,:)' vZ(i,:)'];
            transformCoordAxisRef(i,:)=refInterPolarAxisOrigin+squeeze(recenterCoord(i,1,:))'*refBaseSwitchMatrix;
        end
        
        [azimuth{frameIdx},elevation{frameIdx},rho{frameIdx}]=cart2sph(transformCoordPoleRef(:,1),transformCoordPoleRef(:,2),transformCoordPoleRef(:,3));
        
        % Set each particle in a single cartesian frame of reference again.
%         invertCoord=ones(size(closestPoleIdx));
%         invertCoord(closestPoleIdx==1)=1;
%         movieInfoSpindle(frameIdx).xCoord=[invertCoord.*transformCoord(:,1)+bestPoles(:,1)  pMovieInfo(frameIdx).xCoord(:,2)];
%         movieInfoSpindle(frameIdx).yCoord=[             transformCoord(:,2)+bestPoles(:,2)  pMovieInfo(frameIdx).yCoord(:,2)];
%         movieInfoSpindle(frameIdx).zCoord=[invertCoord.*transformCoord(:,3)+bestPoles(:,3)  pMovieInfo(frameIdx).zCoord(:,2)];

        movieInfoSpindle(frameIdx).xCoord=[transformCoordAxisRef(:,1)/anis(1)  pMovieInfo(frameIdx).xCoord(:,2)];
        movieInfoSpindle(frameIdx).yCoord=[transformCoordAxisRef(:,2)/anis(2)  pMovieInfo(frameIdx).yCoord(:,2)];
        movieInfoSpindle(frameIdx).zCoord=[transformCoordAxisRef(:,3)/anis(3)  pMovieInfo(frameIdx).zCoord(:,2)];
        
    end
end
