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
ip.addParamValue('ROI',[], @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.addParamValue('printAll', false, @islogical);
ip.addParamValue('subDirectory','detections', @ischar);
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

if(p.printAll)
    outputDirDetect=[MD.outputDirectory_ filesep p.subDirectory filesep p.type];
    mkdir(outputDirDetect);
    mkdir([outputDirDetect filesep 'mask']);
end

% find mask offset (WARNING works only for cubic mask)
[maskMinX,maskMinY,maskMinZ]=ind2sub(size(p.ROI), find(p.ROI,1));
[maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(p.ROI), find(p.ROI,1,'last'));

parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(p.channel).loadStack(timePoint));
    volSize=size(vol);
    lab=[];
    if(~isempty(p.ROI))
        tmp=nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
        tmp(:)=vol(p.ROI>0);
        vol=tmp;
    end
    switch p.type
      case 'watershedApplegate'
          filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
            [movieInfo(frameIdx),lab]=detectComets3D(filterVol,p.waterStep,p.waterThresh,[1 1 1]);
      case 'watershedApplegateAuto'
            filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
            thresh=QDApplegateThesh(filterVol,p.showAll);
            [movieInfo(frameIdx),lab]=detectComets3D(filterVol,p.waterStep,thresh,[1 1 1]);
      case 'bandPassWatershed'
        filterVol=filterGauss3D(vol,p.highFreq)-filterGauss3D(vol,p.lowFreq);
        label=watershed(-filterVol); label(filterVol<p.waterThresh)=0;lab=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,filterVol);
      case 'watershedMatlab'
        label=watershed(-vol); label(vol<p.waterThresh)=0;[dummy,nFeats]=bwlabeln(label);
        lab=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,vol);
      case 'markedWatershed'
        lab=markedWatershed(vol,scales,p.waterThresh);
        movieInfo(frameIdx)=labelToMovieInfo(labels{frameIdx},vol);
      case {'pointSource','pointSourceAutoSigma'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
       lab=double(mask); % adjust labellabel
        movieInfo(frameIdx)=labelToMovieInfo(double(mask),vol);
      case {'pointSourceLM','pointSourceAutoSigmaLM'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        lab=double(mask); % adjust label
        movieInfo(frameIdx)=pointCloudToMovieInfo(imgLM,vol);  
      case 'pSAutoSigmaMarkedWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        wat=markedWatershed(vol,scales,0);wat(mask==0)=0;       
       lab=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case 'pSAutoSigmaWatershed'
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        wat=watershed(-vol.*mask);wat(mask==0)=0;
       lab=double(wat); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(wat),vol);
      case {'pointSourceFit','pointSourceAutoSigmaFit'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        lab=double(mask);%.*imgLoG;
        %labels{frameIdx}=double(mask); % adjust label
      case {'pointSourceAutoSigmaMixture'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'FitMixtures', true, varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
       lab=double(mask); % adjust label
      case {'pointSourceAutoSigmaFitSig'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,'Mode','xyzAcsr',varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
       lab=double(mask); % adjust label
      otherwise 
        disp('Unsupported detection method.');
        disp('Supported method:');
        disp('\twatershed');
    end
    
    if(~isempty(p.ROI))
        tmplab=zeros(volSize);
        tmplab(p.ROI>0)=lab;
        lab=tmplab; 
        movieInfo(frameIdx).xCoord(:,1)=movieInfo(frameIdx).xCoord(:,1)+maskMinY-1;
        movieInfo(frameIdx).yCoord(:,1)=movieInfo(frameIdx).yCoord(:,1)+maskMinX-1;
        movieInfo(frameIdx).zCoord(:,1)=movieInfo(frameIdx).zCoord(:,1)+maskMinZ-1;
    end
    
    if p.showAll
        figure()
        vol=double(MD.getChannel(p.channel).loadStack(timePoint));
        imseriesmaskshow(vol,logical(labels{frameIdx}));
%         figure()
%         stackShow(vol,'overlay',labels{frameIdx});
    end
    if(p.printAll)
        stackWrite(uint8(255*lab/max(lab(:))),[outputDirDetect filesep 'mask' filesep 'detect_T_' num2str(timePoint,'%05d') '.tif']);
    end
    labels{frameIdx}=lab;
    
end 

if(p.printAll)
 mkdir([outputDirDetect filesep 'amiraVertex']);
 amiraWriteMovieInfo([outputDirDetect filesep 'amiraVertex' filesep 'detect.am'],movieInfo,'scales',[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_]);
 save([outputDirDetect filesep 'detection.mat'],'movieInfo');
end


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

function mkdir(path)
system(['mkdir -p ' path]);


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