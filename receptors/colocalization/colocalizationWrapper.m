function movieData = colocalizationWrapper(movieData, paramsIn)
%COLOCALIZATIONMOVIE applies a method of colocalization on two channels of an image set
% 
% movieData = colocalizationWrapper(movieData, paramsIn)
%
% Applies a method of colocalization (dependant on how objects are detected
% ) to two channels of an image set and outputs average and sometimes 
% individual colocalization measures
%
% Input:
%   movieData- A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described below.
%       ChannelRef- reference channel. For point2point and continuum2continuum
%       methods, this refers to what fraction of the objects in ChannelObs colocalize
%       with objects in ChannelRef. For point2continuum and point2object, 
%       this is always the point channel (Punctate Channel)
%
%       ChannelObs- observation channel.
%
%       ChannelMask- if masking process is used, indicate here which channel
%       was used for masking
%
%       SearchRadius- distance threshold used for colocalization used in
%       point2point and point2continuum methods
%
%       RandomRuns - number of times positions
%
%       MethodIdx - indicate here which method is to be used. Numbering is
%       as follows: 1-point2point, 2-point2continuum, 3-continuum2continuum
%       4- point2object
%   
% Output: Various colocalization measures dependant on colocalization
% method used for single or stack of images. See specific method function
% for respective outputs. Outputs are saved in unique folder based on
% method (Ex. ColocalizationPt2Cnt/colocalInfo.mat)
%
%   
% Anthony Vega 09/2014

%% Input
% Will need to take previous process outputs from detection and masking 
% processes
%Check that input object is a valid moviedata 
    if nargin < 1 || ~isa(movieData,'MovieData')
        error('The first input argument must be a valid MovieData object!')
    end

    if nargin < 2
        paramsIn = [];
    end
%Get the indices of any previous colocalization processes from this function                                                                              
    iProc = movieData.getProcessIndex('ColocalizationProcess',1,0);

%If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(movieData.processes_)+1;
        movieData.addProcess(ColocalizationProcess(movieData,movieData.outputDirectory_));                                                                                                 
    end


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
    

    %Define which process was masking process
    try
        iM = movieData.getProcessIndex('MaskProcess',Inf,0); 
        inMaskDir = movieData.processes_{iM}.outFilePaths_(p.ChannelMask); 
        maskNames = movieData.processes_{iM}.getOutMaskFileNames(p.ChannelMask);
    catch
        % Try to use imported cell mask if no MaskProcess, Kevin Nguyen 7/2016
        iM = movieData.getProcessIndex('ImportCellMaskProcess',Inf,0); 
        inMaskDir = movieData.processes_{iM}.outFilePaths_{p.ChannelMask}; 
        inMaskDir = fileparts(inMaskDir);
        inMaskDir = {inMaskDir}; % Below requires a cell
        maskNames = {{['cellMask_channel_',num2str(p.ChannelMask),'.tif']}}; 
    end

    
    if p.MethodIndx ~= 3
    %load ref detection data
        try
            iD = movieData.getProcessIndex('ThreshLocMaxProcess',Inf,0);
            inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelRef); 
        catch 

            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelRef); 
        end

        load(inDetectDir{1});
    end
    
    if p.MethodIndx == 1 || p.MethodIndx == 4
        movieInfoRef = movieInfo;
        inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelObs);
        load(inDetectDir{1});
        movieInfoObs = movieInfo;
        pixelIndxObs = pixelIndx;
    end
    

%% Run Colocalization Analysis
switch p.MethodIndx
    case 1 %Point2Point

        [cT,cNull,cCritical,estimatorC,estimatorM]=deal(zeros(length(movieInfoRef),1));
        for iImage = 1:nImages
            detectionRef = movieInfoRef(iImage); 
            detectionObs = movieInfoObs(iImage); 
            currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
            [cT(iImage,:), cNull(iImage,:), cCritical(iImage,:),estimatorM(iImage,:),...
                estimatorC(iImage,:)] = colocalMeasurePt2Pt(detectionRef,detectionObs, p.SearchRadius,currMask);
        end
        mkdir([p.OutputDirectory 'ColocalizationPt2Pt' ])
        save([p.OutputDirectory 'ColocalizationPt2Pt/colocalInfo.mat'],'cT','cNull','cCritical','estimatorC','estimatorM');
    case 2 %Point2Continuum

        [ratioInd,localInd,randRatioInd]=deal(cell(nImages,1));
        start  = movieData.processes_{iD}.funParams_.firstImageNum;
        finish = movieData.processes_{iD}.funParams_.lastImageNum;
        for iImage = start:finish

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
            ratioInd{iImage,:},localInd{iImage,:},randRatioInd{iImage,:},clusterDensity(iImage),cellIntensity(iImage,1)] = colocalMeasurePt2Cnt(p.SearchRadius,...
            p.RandomRuns,detectionData,imageRef,imageObs,currMask);

        end
        mkdir([p.OutputDirectory 'ColocalizationPt2Cnt' ])
        save([p.OutputDirectory 'ColocalizationPt2Cnt/colocalInfo.mat'],'ratioAve','localAve','bgAve','randRatioAve','ratioInd','localInd','randRatioInd','clusterDensity','cellIntensity');
    case 3 %Continuum2Continuum
        for iImage = 1:nImages
            %load image in reference channel
            imageRef = movieData.channels_(p.ChannelRef).loadImage(iImage);
            %load image in observed channel
            imageObs = movieData.channels_(p.ChannelObs).loadImage(iImage);
            %Load the mask for this frame/channel
            currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
            [R(iImage,:),randR(iImage,:)] = colocalMeasureCnt2Cnt(imageRef,imageObs,currMask);
        end
        mkdir([p.OutputDirectory 'ColocalizationCnt2Cnt' ])
        save([p.OutputDirectory 'ColocalizationCnt2Cnt/colocalInfo.mat'], 'R', 'randR');
    case 4 %ColocPt2Obj
        start  = movieData.processes_{iD}.funParams_.firstImageNum;
        finish = movieData.processes_{iD}.funParams_.lastImageNum;
        for iImage = start:finish
            % Load object image
                imageObj = movieData.channels_(p.ChannelObs).loadImage(iImage);        
            %Load detection data
                detectionData = movieInfoRef(iImage); 
            %Load segmentation data
                maskBlobs = zeros(size(imageObj,1),size(imageObj,2));
                for m =1:length(pixelIndxObs{1})
                maskBlobs(pixelIndxObs{1}{m}) = 1;
                end
                bw = bwmorph(maskBlobs,'close');
            %Load the mask for this frame/channel
                currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);

        % Run analysis        
        [colocFract(iImage,:),colocFractRand(iImage,:),colocS{iImage},colocRand{iImage},segProps{iImage}] = colocalMeasurePt2Obj(detectionData,bw,currMask,imageObj);
        end
        mkdir([p.OutputDirectory 'ColocalizationPt2Obj' ])
        save([p.OutputDirectory 'ColocalizationPt2Obj/colocalInfo.mat'],'colocFract','colocFractRand','colocS','colocRand','segProps');
    otherwise
        error('Invalid Method Index!');
end
%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')
end