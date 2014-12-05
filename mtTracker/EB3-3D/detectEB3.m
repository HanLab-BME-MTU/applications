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

if strcmp(ip.Results.type,'pointSourceAutoSigma')||strcmp(ip.Results.type,'pointSourceAutoSigmaFit')
    scales=getGaussianPSFsigmaFrom3DData(double(MD.getChannel(ip.Results.channel).loadStack(1)));
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
        [movieInfo(frameIdx),labels{frameIdx}]=detectComets3D(vol,ip.Results.waterStep,ip.Results.waterThresh,[1 1 1]);
      case 'watershedMatlab'
        label=watershed(-vol); label(vol<ip.Results.waterThresh)=0;[dummy,nFeats]=bwlabeln(label);
        labels{frameIdx}=label;
        movieInfo(frameIdx)=labelToMovieInfo(label,vol);
      case 'markedWatershed'
        [movieInfo(frameIdx),labels{frameIdx}]=markedWatershed(vol,scales,ip.Results.waterThresh);
      case {'pointSource','pointSourceAutoSigma'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        labels{frameIdx}=double(mask); % adjust label
        movieInfo(frameIdx)=labelToMovieInfo(double(mask),vol)
      case {'pointSourceFit','pointSourceAutoSigmaFit'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        movieInfo(frameIdx)= pstructToMovieInfo(pstruct);
        labels{frameIdx}=double(mask); % adjust label
      otherwise 
        disp('Unsupported detection method.');
        disp('Supported method:');
        disp('\twatershed');
    end
    
    if ip.Results.showAll
        figure()
        stackShow(vol,'overlay',labels{frameIdx});
        figure()
        imseriesmaskshow(vol,labels{frameIdx});
    end
end  

function movieInfo= labelToMovieInfo(label,vol)
[feats,nFeats] = bwlabeln(label);
featsProp = regionprops(feats,vol,'Area','WeightedCentroid','MeanIntensity','MaxIntensity','PixelValues');

% centroid coordinates with 0.5 uncertainties
tmp = vertcat(featsProp.WeightedCentroid);
xCoord = [tmp(:,1) 0.5*ones(nFeats,1)]; yCoord = [tmp(:,2) 0.5*ones(nFeats,1)]; zCoord = [tmp(:,3) 0.5*ones(nFeats,1)];
amp=[vertcat(featsProp.MaxIntensity) 0.5*ones(nFeats,1)];

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


