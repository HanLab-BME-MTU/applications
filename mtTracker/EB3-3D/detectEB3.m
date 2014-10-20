function [movieInfo,labels] = detectEB3(MD,varargin)
% SUMMARY: Wrap multiple detection function with standardized input and output for 
%  - easier switching between multiple detector type in a given work flow
%  - systematic comparison
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

if strcmp(ip.Results.type,'pointSourceAutoSigma')
    scales=getGaussianPSFsigmaFrom3DData(double(MD.getChannel(ip.Results.channel).loadStack(1)));
else
    scales=ip.Results.scales;
end 

labels=cell(1,numel(processFrames));
movieInfo=cell(1,numel(processFrames));


parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(ip.Results.channel).loadStack(timePoint));

    switch ip.Results.type
      case 'watershedApplegate'
        [movieInfo,labels{frameIdx}]=detectComets3D(vol,ip.Results.waterStep,ip.Results.waterThresh,[1 1 1]);
      case 'watershedMatlab'
        movieInfo{frameIdx}=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
        label=watershed(-vol); label(vol<ip.Results.waterThresh)=0;[dummy,nFeats]=bwlabeln(label);
        featProp = regionprops(dummy,vol,'WeightedCentroid'); %'Extrema'
        temp = vertcat(featProp.WeightedCentroid);
        movieInfo.yCoord = 0.5*ones(nFeats,3); movieInfo.xCoord = 0.5*ones(nFeats,3); movieInfo.zCoord = 0.5*ones(nFeats,3);
        movieInfo.zCoord(:,1) = temp(:,3); movieInfo.yCoord(:,1) = temp(:,2); movieInfo.xCoord(:,1) = temp(:,1);
        labels{frameIdx}=label;
      case 'markedWatershed'
        [movieInfo,labels{frameIdx}]=markedWatershed(vol,scales,ip.Results.waterThresh);
      case {'pointSource','pointSourceAutoSigma'}
        [pstruct,mask,imgLM,imgLoG]=pointSourceDetection3D(vol,scales,varargin{:});
        labels{frameIdx}=double(mask);
        movieInfo{frameIdx}=pstruct;
      otherwise 
        disp('Unsupported detection method.');
        disp('Supported method:');
        disp('\twatershed');
    end
    
    if ip.Results.showAll
        figure()
        stackShow(vol,'overlay',label);
        figure()
        imseriesmaskshow(vol,label);
    end
end 

% $$$ for frameIndex=ip.Results.processFrames 
% $$$ 
% $$$ if ~strcmp(ip.Results.type,'watershed')
% $$$     [movieInfo,label]=detectComets3D(vol,10,120,[1 1 1]);
% $$$ end {