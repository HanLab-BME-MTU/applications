function [movieInfo,labels] = detectEB3(MD,varargin)
% SUMMARY: Wrap multiple detection function with standardized input and output. Why not using Virtual class ? B/c:
%  - basic and readable switching between multiple detector type in an anylisis process
%  - help systematic comparison
%  - easy detection function inventory
% Philippe Roudot 2014  

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('processFrames',[], @isnumeric);
ip.addParamValue('waterThresh', 120, @isnumeric);
ip.addParamValue('waterStep',10, @isnumeric);
ip.addParamValue('lowFreq',3, @isnumeric);
ip.addParamValue('highFreq',1, @isnumeric);
ip.addParamValue('scales',1.5, @isnumeric);
ip.addParamValue('Alpha',0.05, @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.addParamValue('type', 'watershedApplegate',  @ischar);
ip.parse(MD, varargin{:});

processFrames=[];
if isempty(ip.Results.processFrames)
    processFrames=1:numel(MD.getChannel(ip.Results.channel).getImageFileNames);
else
    processFrames=ip.Results.processFrames;
end

if (strcmp(ip.Results.type,'pointSourceAutoSigma') ...
    ||strcmp(ip.Results.type,'pointSourceAutoSigmaFit') ...
    ||strcmp(ip.Results.type,'pSAutoSigmaMarkedWatershed') ...
    ||strcmp(ip.Results.type,'pointSourceAutoSigmaMixture') ... 
    ||strcmp(ip.Results.type,'pointSourceAutoSigmaLM') ...     
    ||strcmp(ip.Results.type,'pointSourceAutoSigmaFitSig') ... 
    ||strcmp(ip.Results.type,'pSAutoSigmaWatershed'))
    volList=[];
    for i=1:5
        volList=[volList double(MD.getChannel(ip.Results.channel).loadStack(i))];
    end
    scales=getGaussianPSFsigmaFrom3DData(volList);
    disp(['Estimed scales: ' num2str(scales)]);
else
    scales=ip.Results.scales;
end 

labels=cell(1,numel(processFrames));
movieInfo(numel(processFrames),1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);

parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(ip.Results.channel).loadStack(timePoint));
    switch ip.Results.type
      case 'watershedApplegate'
            filterVol=filterGauss3D(vol,ip.Results.highFreq)-filterGauss3D(vol,ip.Results.lowFreq);
            [movieInfo(frameIdx),labels{frameIdx}]=detectComets3D(filterVol,ip.Results.waterStep,ip.Results.waterThresh,[1 1 1]);
      case 'watershedApplegateAuto'
            filterVol=filterGauss3D(vol,ip.Results.highFreq)-filterGauss3D(vol,ip.Results.lowFreq);
            thresh=QDApplegateThesh(filterVol,ip.Results.showAll);
            [movieInfo(frameIdx),labels{frameIdx}]=detectComets3D(filterVol,ip.Results.waterStep,thresh,[1 1 1]);
      case 'bandPassWatershed'
        filterVol=filterGauss3D(vol,ip.Results.highFreq)-filterGauss3D(vol,ip.Results.lowFreq);
        label=watershed(-filterVol); label(filterVol<ip.Results.waterThresh)=0;labels{frameIdx}=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,filterVol);
      case 'watershedMatlab'
        label=watershed(-vol); label(vol<ip.Results.waterThresh)=0;[dummy,nFeats]=bwlabeln(label);
        labels{frameIdx}=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,vol);
      case 'markedWatershed'
        [labels{frameIdx}]=markedWatershed(vol,scales,ip.Results.waterThresh);
        movieInfo(frameIdx)=labelToMovieInfo(labels{frameIdx},vol);
      case {'pointSource','pointSourceAutoSigma'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        labels{frameIdx}=double(mask); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(mask),vol);
      case {'pointSourceLM','pointSourceAutoSigmaLM'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        labels{frameIdx}=double(mask); % adjust label
        movieInfo(frameIdx)=pointCloudToMovieInfo(imgLM,vol);  
      case 'pSAutoSigmaMarkedWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        wat=markedWatershed(vol,scales,0);wat(mask==0)=0;       
        labels{frameIdx}=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case 'pSAutoSigmaWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        wat=watershed(-vol.*mask);wat(mask==0)=0;
        labels{frameIdx}=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case {'pointSourceFit','pointSourceAutoSigmaFit'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      case {'pointSourceAutoSigmaMixture'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'FitMixtures', true, varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      case {'pointSourceAutoSigmaFitSig'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mode','xyzAcsr',varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      otherwise 
        disp('Unsupported detection method.');
        disp('Supported method:');
        disp('\twatershed');
    end
    
    if ip.Results.showAll
        figure()
        imseriesmaskshow(vol,labels{frameIdx});
        figure()
        stackShow(vol,'overlay',labels{frameIdx});
    end
    
end  

function movieInfo= labelToMovieInfo(label,vol)
[feats,nFeats] = bwlabeln(label);
featsProp = regionprops(feats,vol,'Area','WeightedCentroid','MeanIntensity','MaxIntensity','PixelValues');

% centroid coordinates with 0.5 uncertainties
tmp = vertcat(featsProp.WeightedCentroid)-1;
xCoord = [tmp(:,1) 0.5*ones(nFeats,1)]; yCoord = [tmp(:,2) 0.5*ones(nFeats,1)]; zCoord = [tmp(:,3) 0.5*ones(nFeats,1)];
amp=[vertcat(featsProp.MaxIntensity) 0.5*ones(nFeats,1)];

% u-track formating
movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
movieInfo.amp=amp;
movieInfo.int=amp;

function movieInfo= pointCloudToMovieInfo(imgLM,vol)
    lmIdx = find(imgLM~=0);
    [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    N=length(lmy)
    % centroid coordinates with 0.5 uncertainties
    xCoord = [lmx-1 0.5*ones(N,1)]; yCoord = [lmy-1 0.5*ones(N,1)]; zCoord = [lmz-1 0.5*ones(N,1)];
    amp=[vol(lmIdx) 0.5*ones(N,1)];

    % u-track formating
    movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
    movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
    movieInfo.amp=amp;
    movieInfo.int=amp;



function movieInfo= pstructToMovieInfo(pstruct)
movieInfo.xCoord = [pstruct.x'-1 pstruct.x_pstd'];
movieInfo.yCoord = [pstruct.y'-1 pstruct.y_pstd'];
movieInfo.zCoord = [pstruct.z'-1 pstruct.z_pstd'];
movieInfo.amp = [pstruct.A' pstruct.A_pstd'];
movieInfo.int= [pstruct.A' pstruct.A_pstd'];

function threshNoise= QDApplegateThesh(filterDiff,show)
        % Perform maximum filter and mask out significant pixels
        thFilterDiff = locmax3d(filterDiff,1);
        threshold = thresholdOtsu(thFilterDiff)/3 + ...
            thresholdRosin(thFilterDiff)*2/3;
        std=nanstd(filterDiff(thFilterDiff<threshold));
        threshNoise= 3*std;

        if(show)
            figure();hist(thFilterDiff,100),vline([threshNoise, threshold],['-b','-r']);
        end 