function movieData = multiThreshMovie(movieData,paramsIn)
%MULTITHRESHMOVIE applies multi-otsu thresholding to separate the foreground and background of an image
%
% movieData = multiThreshMovie(movieData,paramsIn)
%
% Input:
% movieData - A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described below.
%       ChannelIndex - Number of the channel which will be used for masking
%
%       ThresholdValue - optional input for user chosen threshold value to
%       mask image
%
%       MaxJump - If a threshold is not found for an image in a stack and
%       MaxJump has a value greater than 0, the threshold used in the
%       previous image will be used.
%
%       GaussFilterSigma- the sigma value of the gaussian smoothing filter
%       to be used on the image. Default is no smoothing.
%
%       incThreshold - since the function separates an image into multiple
%       tiers of intensity to find a threshold and the edge of the 
%       foreground maybe poorly defined, it may occur that the
%       threshold chosen is off by a tier. Enter the image number in the
%       stack (or simply 1 if it's a single image) to increase the
%       threshold by one tier.
%
%       decThreshold - similar to incThreshold, enter the image number (or
%       simply 1 if it's a single image) to decrease the chosen threshold
%       by one tier
%
% Output
%       Mask - binary masks for all images are saved in a Mask folder in
%       sub-folders dependant on the channel chosen for masking
%
%       Threshold values- values stored in an array indicating the
%       thresholds used for each image
%
% Anthony Vega 09/2014   
%   
%% ----- Parameters ----- %%

pString = 'mask_'; %Prefix for saving masks to file
pfName = 'threshold_values_for_channel_'; %Prefix for saving threshold values to file. Actual file name will have channel number appended.
dName = 'masks_for_channel_';%String for naming the mask directories for each channel

%% ----------- Input ----------- %%


%Check that input object is a valid moviedata TEMP
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end


%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('MultiThreshProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ThresholdProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

thresProc= movieData.processes_{iProc};


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

nChan = numel(movieData.channels_);

if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end


nChanThresh = length(p.ChannelIndex);
if ~isempty(p.ThresholdValue)
    if length(p.ThresholdValue) == 1
        p.ThresholdValue = repmat(p.ThresholdValue,[1 nChanThresh]);
    elseif length(p.ThresholdValue) ~= nChanThresh
        error('If you specify a threshold value, you must either specify one value to use for all channels, or 1 value per channel!');    
    end
end


%% --------------- Init ---------------%%

disp('Starting thresholding...')

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    if isempty(p.ProcessIndex)
        inFilePaths{1,i} = movieData.getChannelPaths{i};
    else
       inFilePaths{1,i} = movieData.processes_{p.ProcessIndex}.outFilePaths_{1,i}; 
    end
end
thresProc.setInFilePaths(inFilePaths);

%Set up the mask directories as sub-directories of the output directory
for j = 1:nChanThresh;
    
    %Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];    
    %Save this in the process object
    thresProc.setOutMaskPath(p.ChannelIndex(j),currDir);
   
    %Check/create directory
    mkClrDir(currDir)               
end

allThresholds = cell(nChanThresh,1);

imageFileNames = movieData.getImageFileNames(p.ChannelIndex);


nImages = movieData.nFrames_;   
nImTot = nImages * nChanThresh;

%Get mask and image directories
maskDirs  = thresProc.outFilePaths_(p.ChannelIndex);
imDirs  = movieData.getChannelPaths(p.ChannelIndex);
    
%% ----- Thresholding ----- %%

if ~p.BatchMode && feature('ShowFigureWindows')
    wtBar = waitbar(0,['Please wait, thresholding channel ' num2str(p.ChannelIndex(1)) ' ...']);        
else
    wtBar = -1;
end        


for iChan = 1:nChanThresh
        
        
    if ishandle(wtBar)        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, thresholding channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end        
    disp(['Thresholding images for channel # ' num2str(p.ChannelIndex(iChan)) ' : '])
    disp(inFilePaths{1,p.ChannelIndex(iChan)})
    disp('Masks will be stored in directory :')
    disp(maskDirs{iChan})
    
    %Initialize vector for threshold values
    allThresholds{iChan} = nan(nImages,1);
    
    for iImage = 1:nImages
        

        %Load the current image
        if isempty(p.ProcessIndex)
            currImage = movieData.channels_(p.ChannelIndex(iChan)).loadImage(iImage);
        else
            currImage = movieData.processes_{p.ProcessIndex}.loadOutImage(p.ChannelIndex(iChan),iImage);
        end

        %Filter image before thresholding if requested
        if p.GaussFilterSigma > 0
            currImage = filterGauss2D(double(currImage),p.GaussFilterSigma);
        end

        
        if isempty(p.ThresholdValue)
            try
                    fixFrameDown = [];
                    fixFrameUp = [];
                if find(iImage == p.incThreshold)
                    fixFrameUp =1;
                elseif find(iImage == p.decThreshold)
                    fixFrameDown = 1;
                end
                currThresh = calcCellBoundaryImage(currImage,fixFrameUp,fixFrameDown);
                clear fixFrameUp fixFrameDown
            catch %#ok<CTCH>
                %If auto-threshold selection fails, and jump-correction is
                %enabled, force use of previous threshold
                if p.MaxJump > 0
                    currThresh = Inf;
                else
                    if ishandle(wtBar)
                        warndlg(['Could not automatically select a threshold in frame ' ...
                        num2str(iImage) '! Try specifying a threshold level, or enabling the MaxJump option!']);
                        close(wtBar)
                    end                                                
                    error(['Could not automatically select a threshold in frame ' ...
                        num2str(iImage) '! Try specifying a threshold level, or enabling the MaxJump option!']);
                    
                        
                end
            end
        else            
            currThresh = p.ThresholdValue(iChan);
        end
        
        if currThresh == Inf%p.MaxJump > 0
            %Check the threshold
            if iImage == 1
                allThresholds{iChan}(iImage) = currThresh; %Nothing to compare 1st frame to
            else
                if abs(currThresh / allThresholds{iChan}(find(~isnan(allThresholds{iChan}),1,'last'))-1) > (p.MaxJump-1)
                    %If the change was too large, don't store this threshold
                    %and use the most recent good value
                    allThresholds{iChan}(iImage) = NaN;
                    currThresh = allThresholds{iChan}(find(~isnan(allThresholds{iChan}),1,'last'));
                else
                    allThresholds{iChan}(iImage) = currThresh;
                end 
            end
        else
            allThresholds{iChan}(iImage) = currThresh;
        end
        
        %Apply the threshold to create the mask
        imageMask = currImage > currThresh;
        %Test refining images
        test = imfill(imageMask,'holes');
        [L, ~] = bwlabel(test, 8);
        STATS = regionprops(L, 'Area');
        areaInfo = cat(1,STATS.Area);
        sThreshold = find(areaInfo >= (max(areaInfo)/2));
        for i =1:length(sThreshold)
            maskTemp(:,:,i) = (L ==sThreshold(i));
            maskTemp(:,:,1) = maskTemp(:,:,1)+maskTemp(:,:,i);
        end
        imageMask = maskTemp(:,:,1); 

        %write the mask to file                    
        imwrite(imageMask,[maskDirs{iChan} filesep pString imageFileNames{iChan}{iImage}]);
        
        if ishandle(wtBar) && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end
                
    
    end    
   
    
end

if ishandle(wtBar)
    close(wtBar)
end


%% ------ Finish - Save parameters and movieData ----- %%


%Save the threshold values to the analysis directory as seperate files for
%each channel
for i = 1:nChanThresh    
    thresholdValues = allThresholds{i}; %#ok<NASGU>
    save([p.OutputDirectory filesep pfName num2str(p.ChannelIndex(i)) '.mat'],'thresholdValues');
end



thresProc.setDateTime;
movieData.save; %Save the new movieData to disk


disp('Finished thresholding!')

end