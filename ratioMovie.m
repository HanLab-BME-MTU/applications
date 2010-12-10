function movieData = ratioMovie(movieData,paramsIn)
%RATIOMOVIE creates a new ratio channel by dividing one movie channel by another
% 
% movieData = ratioMovie(movieData)
% 
% movieData = ratioMovie(movieData,'OptionName',optionValue)
%
%
% This function divides the images in one channel by those in another and
% stores the resulting images in a new directory.
% 
% 
% Input:
% 
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the corrected images to.
%       Corrected images for different channels will be saved as
%       sub-directories of this directory. If not input, the corrected
%       images will be saved to the same directory as the movieData, in a
%       sub-directory called "bleedthrough_corrected_images"
%
%       ('ChannelIndex'-> Positive integer vector, 1x2) The integer index
%       of the channels to  ratio. This index corresponds to the channel's
%       location in the array movieData.channels_. The first channel
%       specified will be the numerator, and the second will be the
%       denominator.
%       Optional. If not input, the user will be asked to select the
%       channels, unless BatchMode is enabled, in which case an error will
%       be generated.
%       
%       ('ApplyMasks' -> True/False) If true, pixels which are not masked
%       in EITHER channel will be set to zero in the ratio image (i.e. the
%       intersection of the two masks will be used.) Requires that both
%       channels have masks. Default is true.
%
%       ('CreateMasks' -> True/False)
%       If true, the intersection of the masks from both channels will be
%       used to create a new mask for the numerator channel and these masks
%       will over-write the existing numerator masks.
%       Default is false.
%
%       If either CreateMasks or ApplyMasks are true, the following options
%       can also be set:
%
%               ('SegProcessIndex' -> Positive integer scalar or vector)
%               Optional. This specifies SegmentationProcess(s) to use
%               masks from by its index in the array movieData.processes_;
%               If input as a vector, masks will be used from the process
%               specified by the first element, and if not available for a
%               specific channel, then from the next process etc. If not
%               input, and multiple SegmentationProcesses are present, the
%               user will be asked to select one, unless batch mode is
%               enabled in which case there will be an error.  
%
%       ('MaskChannelIndex' -> Positive integer vector, 1x2) The index of
%       the channel of the masks to use for the numerator (1st element) and
%       denominator (2nd element), if masks are applied or created.
%       Default is to use the same channel as the images.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
% Output:
%   
%   movieData - The updated MovieData object, with the parameters and
%   directories for the transformation stored in it as a process object.
%
%   The ratio images will be written to the folder specified by
%   OutputDirectory. They are saved as floating-point, double precision
%   .mat files.
%
% Hunter Elliott, 11/2009
% Revamped 6/2010
%


%%  --------- Parameters ------- %%

pString = 'ratio_'; %The string to prepend before the ratio image directory & channel name
dName = 'ratio_of_'; %String for naming the directories for each ratio channel.


%% ---------- Input -------------%%

%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous ratio processes from this
%function
iProc = movieData.getProcessIndex('RatioProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(RatioProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%----Defaults----%

if isempty(p.ChannelIndex)
    if ~p.BatchMode
        p.ChannelIndex(1) = selectMovieChannels(movieData,0,'Please select the numerator channel:');
        p.ChannelIndex(2) = selectMovieChannels(movieData,0,'Please select the denominator channel:');
    else
        error('In batch mode, the numerator and denominator channels must be specified!')
    end    
end


if length(p.ChannelIndex) ~=2 
    error('You must specify exactly 2 channels for ratioing: a numerator and a denominator!')
end

%Make sure the background subtraction has been performed
iBSProc = find(cellfun(@(x)(isa(x,'BackgroundSubtractionProcess')),movieData.processes_),1);                          
if ~isempty(iBSProc)
    hasBS = cellfun(@(x)(~isempty(x)),movieData.processes_{iBSProc}.outImagePaths_);
else
   error('Background subtraction has not been run! Please run background subtraction prior to ratioing!')   
end

if ~all(hasBS(p.ChannelIndex))
    error('Both channels to be ratioed must be background subtracted prior to ratioing!')
end

nChan = numel(movieData.channels_);

%Check if bleedthrough correction has been run
iBTCProc = find(cellfun(@(x)(isa(x,'BleedthroughCorrectionProcess')),movieData.processes_),1);                          
if ~isempty(iBTCProc)
    %Check which channels have been transformed
    hasBTC = cellfun(@(x)(~isempty(x)),movieData.processes_{iBTCProc}.outImagePaths_);        
else
    hasBTC = false(1,nChan);
end

%Check if transformation has been run
iXFProc = find(cellfun(@(x)(isa(x,'TransformationProcess')),movieData.processes_),1);                          
if ~isempty(iXFProc)
    %Check which channels have been transformed
    hasXF = cellfun(@(x)(~isempty(x)),movieData.processes_{iXFProc}.outImagePaths_);
else
    hasXF = false(1,nChan);
end


if p.ApplyMasks || p.CreateMasks
        
    %Make sure the move has been segmented

    if isempty(p.SegProcessIndex)    
        if p.BatchMode
            %If batch mode, just get all the seg processes
            p.SegProcessIndex = movieData.getProcessIndex('SegmentationProcess',Inf,0);            
        else
            %Otherwise, ask the user 
            p.SegProcessIndex = movieData.getProcessIndex('SegmentationProcess',1,1);
        end
    end

    if isempty(p.SegProcessIndex) 
        error('This function requires that the input movie has already been segmented - no valid SegmentationProcesses were found!')
    end

    nProc = numel(p.SegProcessIndex);
    hasMasks = false(2,nProc);
    
    %Check every specified process for masks
    for j = 1:nProc

        %Make sure the specified process is a SegmentationProcess
        if ~isa(movieData.processes_{p.SegProcessIndex(j)},'SegmentationProcess')
            error(['The process specified by SegProcessIndex(' num2str(j) ') is not a SegmentationProcess!'])
        end

        %Check which channels have masks from this process
        hasMasks(:,j) = movieData.processes_{p.SegProcessIndex(j)}.checkChannelOutput(p.MaskChannelIndex);        

    end

    %Make sure all the selected channels have foreground masks.
    if any(~sum(hasMasks,2))
        warning('biosensors:backgroundMasks:noFGmasks',...
            'Cannot create / apply masks because some channels do not have masks! Please segment these channels before creating / applying ratio masks!')
    end           
        
    %Get the first seg process with masks for this channel
      %Get the first seg process with masks for this channel
    iP = p.SegProcessIndex(find(hasMasks(1,:),1));
    
    numMaskDir = movieData.processes_{iP}.outMaskPaths_{p.MaskChannelIndex(1)};
    numMaskNames = movieData.processes_{iP}.getOutMaskFileNames(p.MaskChannelIndex(1));        
    
    iP = p.SegProcessIndex(find(hasMasks(2,:),1));
    
    denomMaskDir = movieData.processes_{iP}.outMaskPaths_{p.MaskChannelIndex(2)};        
    denomMaskNames = movieData.processes_{iP}.getOutMaskFileNames(p.MaskChannelIndex(2));
                          
end

%Save these selected channels / parameters in the movieData
movieData.processes_{iProc}.setPara(p)    


%Set up input/output directories
for j = 1:2
    if hasXF(p.ChannelIndex(j))
        %If available, use transformed images
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iXFProc}.outImagePaths_{p.ChannelIndex(j)});
    elseif hasBTC(p.ChannelIndex(j))
        %... or bleedthrough corrected images
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iBTCProc}.outImagePaths_{p.ChannelIndex(j)});
    else
        %Otherwise, use background subtracted.
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iBSProc}.outImagePaths_{p.ChannelIndex(j)});
    end                    

end


%% ------------- Init ------------ %%

disp('Initializing ratioing...')

%Ratios are stored in numerator channel
outDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(1)) '_to_' ...
    num2str(p.ChannelIndex(2))];

%Check/set up output directory
mkClrDir(outDir);

movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex(1),outDir);    

nImages = movieData.nFrames_;

numDir = movieData.processes_{iProc}.inImagePaths_{p.ChannelIndex(1)};
numImNames = movieData.processes_{iProc}.getInImageFileNames(p.ChannelIndex(1));
denomDir = movieData.processes_{iProc}.inImagePaths_{p.ChannelIndex(2)};
denomImNames = movieData.processes_{iProc}.getInImageFileNames(p.ChannelIndex(2));

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nImages))+1) '.f'];


ratMaskDir = [p.OutputDirectory filesep 'ratio_masks'];
mkClrDir(ratMaskDir);

%% ------ Ratio -----%%
% Ratios the channels and writes the resulting ratio images to file


disp('Starting ratioing...')
disp(['Creating ratio images by dividing channel ' numDir ' by channel ' denomDir])
disp(['Resulting images will be written to channel ' outDir])
  
if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, ratioing channel ' ...
        num2str(p.ChannelIndex(1)) ' to channel ' num2str(p.ChannelIndex(2)) '...']);        
end        

for iImage = 1:nImages
    
     
    currNum = imread([numDir filesep numImNames{1}{iImage}]);   
    currDenom = imread([denomDir filesep denomImNames{1}{iImage}]);
    
    %No big deal....
    currRatio = double(currNum) ./ double(currDenom);
    
    %We create our own zero-padded numbering to cover any file-naming
    %problems (e.g. non zero-padded metamorph images!)
    numStr = num2str(iImage,fString);
    
    
   if p.ApplyMasks || p.CreateMasks  %If masks are to be applied, don't include the masked values      
        
        if hasMasks(1)
            currNumMask = imread([numMaskDir filesep numMaskNames{1}{iImage}]);
            if p.ApplyMasks 
                currRatio(~currNumMask(:)) = 0;        
            end
        elseif p.CreateMasks
            currNumMask = true(size(currRatio));
        end
        
        if hasMasks(2)
            currDenomMask = imread([denomMaskDir filesep denomMaskNames{1}{iImage}]);
            if p.ApplyMasks 
                currRatio(~currDenomMask(:)) = 0;        
            end
        elseif p.CreateMasks
            currDenomMask = true(size(currRatio));
        end
        if p.CreateMasks              
            imwrite(currDenomMask & currNumMask,[ratMaskDir filesep 'ratio_mask_' numStr '.tif'])
        end
        
    end    
    
    %Remove any infinities from division-by-zero (this shouldn't happen if
    %the masks are applied and images are good, but let's be realistic here .... )
    currRatio(~isfinite(currRatio(:))) = 0; %#ok<NASGU>
  
    %Save the ratio in double-precision floating point to avoid rounding
    %error
    
  
    save([outDir filesep pString numImNames{1}{iImage}(1:end-4) ...
        '_to_' denomImNames{1}{iImage}(1:end-4) '_' numStr ...
        '.mat'],'currRatio')
    
    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iImage/nImages,wtBar)
    end                            
    
end


if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end



%% ------ Output / Finalization ---- %%

disp('Saving results...')

%Log the correction in the movieData object and save it

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData;

disp('Finished Ratioing!')


