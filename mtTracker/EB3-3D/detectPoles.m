function [poleMovieInfo] = detectPoles(MD,varargin)
% Philippe Roudot 2014  

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('processFrames',[], @isnumeric);
ip.addParamValue('scales',[10 10 10], @isnumeric);
ip.addParamValue('Alpha',0.05, @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.addParamValue('printAll', false, @islogical);
ip.addParamValue('type', 'simplex',  @ischar);
ip.parse(MD, varargin{:});

processFrames=[];
if isempty(ip.Results.processFrames)
    processFrames=1:numel(MD.getChannel(ip.Results.channel).getImageFileNames);
else
    processFrames=ip.Results.processFrames;
end

scales=ip.Results.scales;
poleMovieInfo(numel(processFrames),1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);
movieInfo(numel(processFrames),1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);
parfor frameIdx=1:numel(processFrames)
    timePoint=processFrames(frameIdx);
    disp(['Processing time point ' num2str(timePoint,'%04.f')])
    vol=double(MD.getChannel(ip.Results.channel).loadStack(timePoint));
    ws = ceil(2*scales);

    gx = exp(-(0:ws(1)).^2/(2*scales(1)^2));
    gz = exp(-(0:ws(3)).^2/(2*scales(3)^2));
    fg = conv3fast(vol, gx, gx, gz);
    
    lm=locmaxnd(fg,ceil(scales));
%     lm(1:ws(1),:,:)=0;
%     lm(:,1:ws(2),:)=0;
%     lm(:,:,1:ws(3))=0;
    perc=100;
    notEnoughPoles=true;
    while( notEnoughPoles)
        perc=perc-5;
        percentile=prctile(lm(lm>0),perc);
        notEnoughPoles=sum(sum(sum(lm>=percentile)))<2;
    end
    lm(lm<percentile)=0;
    movieInfo(frameIdx)=pointCloudToMovieInfo(lm,vol);
end  

%% load track results and save them to Amira
if(ip.Results.printAll)
outputDirDetect=[MD.outputDirectory_ filesep 'poles' filesep ip.Results.type '_scale_' ... 
    num2str(scales(1),'%03d') filesep 'poleCandidates'];
mkdir([outputDirDetect filesep 'AmiraPoles']);
amiraWriteMovieInfo([outputDirDetect filesep filesep 'polesCandidates.am'],movieInfo,'scales',ip.Results.scales);
end

%% Track each candidate to filter theire intensit
[gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=candidatePolesTrackingParam();
outputDirTrack=[MD.outputDirectory_ filesep 'poles' filesep ip.Results.type '_scale_' ... 
    num2str(scales(1),'%03d') filesep 'tracks'];

saveResults.dir =  outputDirTrack; %directory where to save input and output
saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved
saveResults=[];
[tracksFinal,kalmanInfoLink,errFlag] = ...
    trackCloseGapsKalmanSparse(movieInfo, ...
    costMatrices,gapCloseParam,kalmanFunctions,...
    probDim,saveResults,verbose);

%% Convert tracks final in a user-friendlier format
tracks=TracksHandle(tracksFinal);

if(ip.Results.printAll)
    mkdir(outputDirTrack);
    save([outputDirTrack filesep 'trackNewFormat.mat'],'tracks')
end 

%% Retrieve innovation matrix 
trackNoiseVar=arrayfun(@(x) kalmanInfoLink(tracks(x).segmentEndFrame).noiseVar(1,1,tracks(x).tracksFeatIndxCG(end)),1:length(tracks))';

%% Save track results to Amira
if(ip.Results.printAll)
    mkdir([outputDirTrack filesep 'AmiraTrack']);
    amiraWriteTracks([outputDirTrack filesep 'AmiraTrack' filesep 'test.am'],tracks,'scales',[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_],'edgeProp',{{'noiseVar',trackNoiseVar}});
end

%% For each frame, select the higher responses,
tracksScore=[tracks.lifetime].*arrayfun(@(x) median(x.A),tracks)';

%% Compute the distance between each candidate (looking for stationary distance maybe ?) 
% tracksMeanPos=[arrayfun(@(x) median(x.x),tracks) arrayfun(@(x) median(x.y),tracks) arrayfun(@(x) median(x.z),tracks)]
% distMatrix=createSparseDistanceMatrix(tracksMeanPos,tracksMeanPos,1000000);
% trackMaxDist=full(max(distMatrix));
% tracksScore=trackMaxDist;


for fIdx=1:numel(processFrames)
    timePoint=processFrames(fIdx);
    tracksOn=([tracks.endFrame]>=timePoint)&(timePoint>=[tracks.startFrame]);
    tracksLocal=tracks(tracksOn);
    relIdx=timePoint-[tracksLocal.startFrame]+1; 
    [~,idx]=sort(tracksScore(tracksOn));
    selectedIdx=[];
    if(sum(tracksOn)>0)
        IdxPole1=(tracksLocal(idx(end)).tracksFeatIndxCG(relIdx(idx(end))));
        selectedIdx=[selectedIdx IdxPole1];
    end;
    if(sum(tracksOn)>1)
        IdxPole2=(tracksLocal(idx(end-1)).tracksFeatIndxCG(relIdx(idx(end-1))));
        selectedIdx=[selectedIdx IdxPole2];
    end;
    MI=movieInfo(timePoint);
    fn=fieldnames(MI);
    for i=1:length(fn) poleMovieInfo(fIdx).(fn{i})=MI.(fn{i})(selectedIdx,:); end;
end  


function movieInfo= labelToMovieInfo(label,vol)
[feats,nFeats] = bwlabeln(label);
featsProp = regionprops(feats,vol,'Area','WeightedCentroiccd','MeanIntensity','MaxIntensity','PixelValues');

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