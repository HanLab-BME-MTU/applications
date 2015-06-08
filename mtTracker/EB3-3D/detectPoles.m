function [movieInfo] = detectPoles(MD,varargin)
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
ip.addParamValue('scales',[10 10 10], @isnumeric);
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


scales=ip.Results.scales;
 

labels=cell(1,numel(processFrames));
movieInfo(numel(processFrames),1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);

parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(ip.Results.channel).loadStack(timePoint));
    ws = ceil(2*scales);

    gx = exp(-(0:ws(1)).^2/(2*scales(1)^2));
    gz = exp(-(0:ws(3)).^2/(2*scales(3)^2));
    fg = conv3fast(vol, gx, gx, gz);
    
    lm=locmaxnd(fg,scales);
    sortedMax=sort((lm(lm>0)),1,'descend');
    lm((lm~=sortedMax(1))&(lm~=sortedMax(2)))=0;
    movieInfo(frameIdx)=pointCloudToMovieInfo(lm,vol);
    
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
    N=length(lmy);
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