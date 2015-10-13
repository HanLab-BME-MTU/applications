function [movieInfo,labels,scales] = detectEB3(MD,varargin)
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
ip.addParamValue('scales',[], @isnumeric);
ip.addParamValue('Alpha',0.05, @isnumeric);
ip.addParamValue('mask',[], @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.addParamValue('printAll', false, @islogical);
ip.addParamValue('subDirectory','EB3', @ischar);
ip.addParamValue('type', 'watershedApplegate',  @ischar);
ip.parse(MD, varargin{:});

p=ip.Results;

processFrames=[];
if isempty(p.processFrames)
    processFrames=1:numel(MD.getChannel(p.channel).getImageFileNames);
else
    processFrames=p.processFrames;
end

if (strcmp(p.type,'pointSourceAutoSigma') ...
    ||strcmp(p.type,'pointSourceAutoSigmaFit') ...
    ||strcmp(p.type,'pSAutoSigmaMarkedWatershed') ...
    ||strcmp(p.type,'pointSourceAutoSigmaMixture') ... 
    ||strcmp(p.type,'pointSourceAutoSigmaLM') ...     
    ||strcmp(p.type,'pointSourceAutoSigmaFitSig') ... 
    ||strcmp(p.type,'pSAutoSigmaWatershed')) ... 
    &&(isempty(p.scales))
    volList=[];
    for i=1:5
        volList=[volList double(MD.getChannel(p.channel).loadStack(i))];
    end
    scales=getGaussianPSFsigmaFrom3DData(volList,'Display',p.showAll);
    disp(['Estimed scales: ' num2str(scales)]);
else
    scales=p.scales;
end 

labels=cell(1,numel(processFrames));
movieInfo(numel(processFrames),1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);

outputDirDetect=[MD.outputDirectory_ filesep p.subDirectory filesep p.type];
mkdir(outputDirDetect);
mkdir([outputDirDetect filesep 'mask']);
parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(p.channel).loadStack(timePoint));
    switch p.type
      case 'watershedApplegate'
          filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
            [movieInfo(frameIdx),labels{frameIdx}]=detectComets3D(filterVol,p.waterStep,p.waterThresh,[1 1 1]);
      case 'watershedApplegateAuto'
            filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
            thresh=QDApplegateThesh(filterVol,p.showAll);
            [movieInfo(frameIdx),labels{frameIdx}]=detectComets3D(filterVol,p.waterStep,thresh,[1 1 1]);
      case 'bandPassWatershed'
        filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
        label=watershed(-filterVol); label(filterVol<p.waterThresh)=0;labels{frameIdx}=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,filterVol);
      case 'watershedMatlab'
        label=watershed(-vol); label(vol<p.waterThresh)=0;[dummy,nFeats]=bwlabeln(label);
        labels{frameIdx}=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,vol);
      case 'markedWatershed'
        [labels{frameIdx}]=markedWatershed(vol,scales,p.waterThresh);
        movieInfo(frameIdx)=labelToMovieInfo(labels{frameIdx},vol);
      case {'pointSource','pointSourceAutoSigma'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,varargin{:});
        labels{frameIdx}=double(mask); % adjust labellabel
        movieInfo(frameIdx)=labelToMovieInfo(double(mask),vol);
      case {'pointSourceLM','pointSourceAutoSigmaLM'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,varargin{:});
        labels{frameIdx}=double(mask); % adjust label
        movieInfo(frameIdx)=pointCloudToMovieInfo(imgLM,vol);  
      case 'pSAutoSigmaMarkedWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,varargin{:});
        wat=markedWatershed(vol,scales,0);wat(mask==0)=0;       
        labels{frameIdx}=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case 'pSAutoSigmaWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,varargin{:});
        wat=watershed(-vol.*mask);wat(mask==0)=0;
        labels{frameIdx}=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case {'pointSourceFit','pointSourceAutoSigmaFit'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
          labels{frameIdx}=double(mask);
        %labels{frameIdx}=double(mask); % adjust label
      case {'pointSourceAutoSigmaMixture'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,'FitMixtures', true, varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      case {'pointSourceAutoSigmaFitSig'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mask',p.mask,'Mode','xyzAcsr',varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      otherwise 
        disp('Unsupported detection method.');
        disp('Supported method:');
        disp('\twatershed');
    end
    
    if p.showAll
        figure()
        imseriesmaskshow(vol,logical(labels{frameIdx}));
%         figure()
%         stackShow(vol,'overlay',labels{frameIdx});
    end
    if(p.printAll)
        stackWrite(labels{frameIdx},[outputDirDetect filesep 'mask' filesep 'detect_T_' num2str(frameIdx,'%05d') '.tif']);
    end
    
end 

if(p.printAll)
 mkdir([outputDirDetect filesep 'amiraVertex']);
 amiraWriteMovieInfo([outputDirDetect filesep 'amiraVertex' filesep 'detect.am'],movieInfo,'scales',[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_]);
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
    N=length(lmy);
    % centroid coordinates with 0.5 uncertainties
    xCoord = [lmx 0.5*ones(N,1)]; yCoord = [lmy 0.5*ones(N,1)]; zCoord = [lmz 0.5*ones(N,1)];
    amp=[vol(lmIdx) 0.5*ones(N,1)];

    % u-track formating
    movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
    movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
    movieInfo.amp=amp;
    movieInfo.int=amp;



function movieInfo= pstructToMovieInfo(pstruct)
movieInfo.xCoord = [pstruct.x' pstruct.x_pstd'];
movieInfo.yCoord = [pstruct.y' pstruct.y_pstd'];
movieInfo.zCoord = [pstruct.z' pstruct.z_pstd'];
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