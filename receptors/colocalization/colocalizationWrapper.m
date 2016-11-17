function movieData = colocalizationWrapper(movieData, paramsIn)
%COLOCALIZATIONWRAPPER applies colocalization method described in colocalMeasurePt2Cnt to a single image or an image set
% 
% movieData = colocalizationWrapper(movieData, paramsIn)
%
% Applies colocalization method described in colocalMeasurePt2Cnt to two 
% channels (punctate and continuum) of an image set and outputs average and
% sometimes individual colocalization measures
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
%   
% Output: See core function, colocalMeasurePt2Cnt
% for specific outputs. Outputs are saved in unique folder 
% ColocalizationPt2Cnt/colocalInfoMN.mat where M = punctate channel and N=
% continuum channel number used in analysis.
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

    %load ref detection data
        try
            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelRef); 
        catch 
            %Insert some error about no detection file
    
        end

        load(inDetectDir{1});

%% Run Colocalization Analysis

        [enrichInd,localInd,randEnrichInd]=deal(cell(nImages,1));
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

            [enrichAve(iImage,:),localAve(iImage,:),bgAve(iImage,:),randEnrichAve(iImage,:),...
            enrichInd{iImage,:},localInd{iImage,:},randEnrichInd{iImage,:},clusterDensity(iImage),cellIntensity(iImage,1)] = colocalMeasurePt2Cnt(p.SearchRadius,...
            p.RandomRuns,detectionData,imageRef,imageObs,currMask);

        end
        mkdir([p.OutputDirectory 'ColocalizationPt2Cnt' ])
        save([p.OutputDirectory 'ColocalizationPt2Cnt/colocalInfo' num2str(p.ChannelRef) num2str(p.ChannelObs) '.mat'],'enrichAve','localAve','bgAve','randEnrichAve','enrichInd','localInd','randEnrichInd','clusterDensity','cellIntensity');


%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')
end