function movieData = colocalizationWrapper(movieData, paramsIn)
%COLOCALIZATIONMOVIE applies a method of colocalization on two channels of an image set
% 
% movieData = colocalizationWrapper(movieData, paramsIn)
%
% Applies a method of colocalization (pt2pt, pt2cnt) to two channels of an
% image set and outputs average and sometimes individual colocalization
% measures
%
% Input:
% Output:
% Anthony Vega 09/2014

%% Input
% Will need to take previous process outputs from detection and masking 
% processes
%Check that input object is a valid moviedata TEMP
    if nargin < 1 || ~isa(movieData,'MovieData')
        error('The first input argument must be a valid MovieData object!')
    end

    if nargin < 2
        paramsIn = [];
    end
%Get the indices of any previous threshold processes from this function                                                                              
    iProc = movieData.getProcessIndex('ColocalizationProcess',1,0);

%If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(movieData.processes_)+1;
        movieData.addProcess(ColocalizationProcess(movieData,movieData.outputDirectory_));                                                                                                 
    end

    colocProc= movieData.processes_{iProc};


%Parse input, store in parameter structure
    p = parseProcessParams(movieData.processes_{iProc},paramsIn);

    nChan = numel(movieData.channels_);

    if max(p.ChannelRef) > nChan || min(p.ChannelRef)<1 || ~isequal(round(p.ChannelRef),p.ChannelRef)
        error('Invalid channel numbers specified! Check ChannelRef input!!')
    end

    if max(p.ChannelObs) > nChan || min(p.ChannelObs)<1 || ~isequal(round(p.ChannelObs),p.ChannelObs)
        error('Invalid channel numbers specified! Check ChannelObs input!!')
    end


% After parameters are chosen, take appropriate channels from movieData
    nImages = movieData.nFrames_;
    colocMethod = colocProc.getMethods(p.MethodIndx).func; %Not sure if necessary
  
    %load masking data
    nChanThresh = length(p.ChannelMask);
    inMaskDir = cell(1,nChanThresh);
    maskNames = cell(1,nChanThresh);
    
    %Define which process was masking process?
    iM = movieData.getProcessIndex('MaskProcess',Inf,0);
    inMaskDir = movieData.processes_{iM}.outFilePaths_(p.ChannelMask); 
    maskNames = movieData.processes_{iM}.getOutMaskFileNames(p.ChannelMask);
    
    if p.MethodIndx ~= 3
    %load ref detection data
    iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
    inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelRef); 
    load(inDetectDir{1});
    end
    if p.MethodIndx == 1
        movieInfoRef = movieInfo;
        inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelObs);
        load(inDetectDir{1});
        movieInfoObs = movieInfo;
    end
    %Create Structure

%% Run Colocalization Analysis
switch p.MethodIndx
    case 1
% %         colocalInfoAve = repmat(struct('ratioAve',[],'localAve',[],'bgAve',[],'randRatioAve',[]),nImages,1);
% %         colocalInfoInd = repmat(struct('ratioInd',[],'localInd',[],'randRatioInd',[]),nImages,1);
        [cT,cNull,cCritical,estimatorC,estimatorM]=deal(zeros(length(movieInfoRef),1));
        for iImage = 1:nImages
            detectionRef = movieInfoRef(iImage); 
            detectionObs = movieInfoObs(iImage); 
            currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
            [cT(iImage,:), cNull(iImage,:), cCritical(iImage,:),estimatorM(iImage,:),...
                estimatorC(iImage,:)] = colocalMeasurePt2Pt(detectionRef,detectionObs, p.SearchRadius,currMask);
        end    
        save(p.OutputDirectory,'cT','cNull','cCritical','estimatorC','estimatorM');
    case 2
%         colocalInfo = repmat(struct('ratioAve',[],'localAve',[],'bgAve',[],'randRatioAve',[],'ratioInd',[],'localInd',[],'randRatioInd',[]),nImages,1);
%         colocalInfoAve = repmat(struct('ratioAve',[],'localAve',[],'bgAve',[],'randRatioAve',[]),nImages,1);
%         colocalInfoInd = repmat(struct('ratioInd',[],'localInd',[],'randRatioInd',[]),nImages,1);
%         []=deal(zeros(length(imageAInfo),1));
        [ratioInd,localInd,randRatioInd]=deal(cell(nImages,1));
        for iImage = 1:nImages
    if ~isempty(find(iImage ==[28] ))
% 
        continue
    end
            %load image in reference channel
            imageRef = movieData.channels_(p.ChannelRef).loadImage(iImage);
            %load image in observed channel
            imageObs = movieData.channels_(p.ChannelObs).loadImage(iImage);

            %Load detection data
            detectionData = movieInfo(iImage);         
            %Load the mask for this frame/channel
            currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
        
            %Run Function

            [ratioAve(iImage,:),localAve(iImage,:),bgAve(iImage,:),randRatioAve(iImage,:),...
            ratioInd{iImage,:},localInd{iImage,:},randRatioInd{iImage,:},clusterDensity(iImage),cellIntensity(iImage,1),backgroundStd(iImage,1)] = colocalMeasurePt2Cnt(p.SearchRadius,...
            p.RandomRuns,detectionData,imageRef,imageObs,currMask);
%             colocalInfoAve(iImage) = colocalFeaturesAve;
%             colocalInfoInd(iImage) = colocalFeaturesInd;
        end
        save(p.OutputDirectory,'ratioAve','localAve','bgAve','randRatioAve','ratioInd','localInd','randRatioInd','clusterDensity','cellIntensity','backgroundStd');
    case 3
        for iImage = 1:nImages
            %load image in reference channel
            imageRef = movieData.channels_(p.ChannelRef).loadImage(iImage);
            %load image in observed channel
            imageObs = movieData.channels_(p.ChannelObs).loadImage(iImage);
            %Load the mask for this frame/channel
            currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
            [R(iImage,:),randR(iImage,:)] = colocalMeasureCnt2Cnt(imageRef,imageObs,currMask);
        end
        save(p.OutputDirectory, 'R', 'randR');
    otherwise
        error('Invalid Method Index!');
end
%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')
end