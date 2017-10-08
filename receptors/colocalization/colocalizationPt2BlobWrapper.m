function movieData = colocalizationPt2BlobWrapper(movieData, paramsIn)
%COLOCALIZATIONPT2BLOBWRAPPER applies colocalization method described in colocalMeasurePt2Cnt to a single image or an image set
% 
% movieData = colocalizationPt2BlobWrapper(movieData, paramsIn)
%
% Applies colocalization method described in colocalMeasurePt2Blob to two 
% channels (punctate and blob(segmentable objects)) of an image set and outputs average and
% sometimes individual colocalization measures
%
% Input:
%   movieData- A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described below.
%       ChannelPt- reference channel. For point2point and continuum2continuum
%       methods, this refers to what fraction of the objects in ChannelObs colocalize
%       with objects in ChannelRef. For point2continuum and point2object, 
%       this is always the point channel (Punctate Channel)
%
%       ChannelBlob- observation channel.
%
%       ChannelMask- if masking process is used, indicate here which channel
%       was used for masking
%
%       SearchRadius- distance threshold used for colocalization used in
%       point2point, point2blob, and point2continuum methods
%
%
%   
% Output: See core function, colocalMeasurePt2Blob
% for specific outputs. Outputs are saved in unique folder 
% ColocalizationPt2Blob/colocalInfoMN.mat where M = punctate channel and N=
% continuum channel number used in analysis.
%
%   
% Anthony Vega 09/2017

%% Input
% Will need to take previous process outputs from detection,segmentation and masking 
% processes
%Check that input object is a valid moviedata 
    if nargin < 1 || ~isa(movieData,'MovieData')
        error('The first input argument must be a valid MovieData object!')
    end

    if nargin < 2
        paramsIn = [];
    end
%Get the indices of any previous colocalization processes from this function                                                                              
    iProc = movieData.getProcessIndex('ColocalizationPt2BlobProcess',1,0);

%If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(movieData.processes_)+1;
        movieData.addProcess(ColocalizationProcess(movieData,movieData.outputDirectory_));                                                                                                 
    end


%Parse input, store in parameter structure
    p = parseProcessParams(movieData.processes_{iProc},paramsIn);

    nChan = numel(movieData.channels_);

    if max(p.ChannelPt) > nChan || min(p.ChannelPt)<1 || ~isequal(round(p.ChannelPt),p.ChannelPt)
        error('Invalid channel numbers specified! Check ChannelPt input!!')
    end

    if max(p.ChannelBlob) > nChan || min(p.ChannelBlob)<1 || ~isequal(round(p.ChannelBlob),p.ChannelBlob)
        error('Invalid channel numbers specified! Check ChannelBlob input!!')
    end


% After parameters are chosen, take appropriate channels from movieData
    nImages = movieData.nFrames_;
    

    %Define which process was masking process
    if ~isempty(p.ChannelMask)
        try
            warning('off', 'lccb:process')
            iM = movieData.getProcessIndex('MaskProcess',1,0); 
            inMaskDir = movieData.processes_{iM}.outFilePaths_(p.ChannelMask); 
            maskNames = movieData.processes_{iM}.getOutMaskFileNames(p.ChannelMask);
            warning('on', 'lccb:process')
        catch
            try
                warning('off', 'lccb:process')
                iM = movieData.getProcessIndex('MultiThreshProcess',1,0); 
                inMaskDir = movieData.processes_{iM}.outFilePaths_(p.ChannelMask); 
                maskNames = movieData.processes_{iM}.getOutMaskFileNames(p.ChannelMask);
                warning('on', 'lccb:process')
            catch
                % Try to use imported cell mask if no MaskProcess, Kevin Nguyen 7/2016
                iM = movieData.getProcessIndex('ImportCellMaskProcess',Inf,0); 
                inMaskDir = movieData.processes_{iM}.outFilePaths_{p.ChannelMask}; 
                inMaskDir = fileparts(inMaskDir);
                inMaskDir = {inMaskDir}; % Below requires a cell
                maskNames = {{['cellMask_channel_',num2str(p.ChannelMask),'.tif']}};
            end
        end
    end
    %load pt channel detection data
        try
            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            inDetectDir = movieData.processes_{iD}.outFilePaths_(p.ChannelPt); 
        catch 
            %Insert some error about no detection file
            error('No detection data found! Detection process must be run first!')
        end
    %Load blob channel segmentation data
        try
            iS = movieData.getProcessIndex('SegmentBlobsPlusLocMaxSimpleProcess',Inf,0);
            inBlobSegDir = movieData.processes_{iS}.outFilePaths_(p.ChannelBlob); 
        catch 
            %Insert some error about no segmentation file
            error('No segmentation data found! Detection process must be run first!')
        end
        load(inDetectDir{1});
        load(inBlobSegDir{1});

%% Run Colocalization Analysis

        [cT,cNull,cCritical]=deal(zeros(nImages,1));
        start  = movieData.processes_{iD}.funParams_.firstImageNum;
        finish = movieData.processes_{iD}.funParams_.lastImageNum;
        for iImage = start:finish
          

            %load image in blob channel
            imageBlob = movieData.channels_(p.ChannelBlob).loadImage(iImage);

            %Load detection data
            detectionData = movieInfo(iImage);         
            %Load the mask for this frame/channel
            if ~isempty(p.ChannelMask)
                currMask = imread([inMaskDir{1} filesep maskNames{1}{iImage}]);
            else
                currMask = [];
            end
            %Get segmenation coordinates from segmentation data 
            [segmentationData] = segmentationCoordinates(imageBlob,maskBlobs,[],[],[],[]);
            
            %Run Function

            [cT(iImage), cNull(iImage), cCritical(iImage)] = colocalMeasurePt2Blob(segmentationData,detectionData, p.SearchRadius, currMask);          

        end
        mkdir([p.OutputDirectory 'ColocalizationPt2Blob' ])
        save([p.OutputDirectory 'ColocalizationPt2Blob/colocalInfo' num2str(p.ChannelPt) num2str(p.ChannelBlob) '.mat'],'cT','cNull','cCritical');


%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')