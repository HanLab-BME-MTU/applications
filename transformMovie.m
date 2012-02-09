function movieData = transformMovie(movieData,paramsIn)

% movieData = transformMovie(movieData)
% 
% movieData = transformMovie(movieData,paramsIn)
% 
% This function performs a spatial transformation on the selected channels
% of the input movie and writes the transformed images to a new channel in
% the movie. The transformation should be saved as a .mat file.
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
%       ('ChannelIndex'-> Positive integer scalar) The integer index of the
%       channel to perform spatial transformation on. This index
%       corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be
%       transformed.
%       
%       ('TransformFilePaths' -> Cell array of Character strings) A cell
%       array specifying The FULL path and filename of the .mat file
%       containing the transform to apply to the images in each channel.
%       Should contain one element for each channel to be transformed. The
%       transform should be of the format used by imtransform.m. If not
%       input, the user will be asked to locate a file containing a
%       transform for each channel, UNLESS batchmode is enabled, in which
%       case an error will be generated.
%
%       ('TransformMasks' -> True/False)
%       If true, the masks for a given channel will also be transformed,
%       and saved to a mask directory for the output channel(s). If true,
%       the specified channels MUST have masks. Default is true.
%   
%       If masks are transformed, these additional options apply:
%       
%               ('SegProcessIndex' -> Positive integer scalar or vector)
%               Optional. This specifies MaskProcess(s) to use
%               masks from by its index in the array movieData.processes_;
%               If input as a vector, masks will be used from the process
%               specified by the first element, and if not available for a
%               specific channel, then from the next process etc. If not
%               input, and multiple MaskProcesses are present, the
%               user will be asked to select one, unless batch mode is
%               enabled in which case there will be an error.  
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed. Default is false.
%
% 
% Output:
%   
%   movieData - The updated MovieData object, with the parameters and
%   directories for the transformation stored in it as a process object.
%
%   The transformed images will be written to the folder specified by
%   OutputDirectory.
%
% Hunter Elliott
% 11/2009
% Revamped 5/2010
%

%% ------ Parameters ----%%

pString = 'xf_'; %The string to prepend before the transformed image directory & channel name
dName = 'transformed_images_for_channel_';%String for naming the directories for each corrected channel

%% ------- Input ------- %%

%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous bleedthrough correction processes from this
%function
iProc = find(cellfun(@(x)(isa(x,'TransformationProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(TransformationProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has been background-subtracted
iBSProc = find(cellfun(@(x)(isa(x,'BackgroundSubtractionProcess')),movieData.processes_),1);                          
if isempty(iBSProc)
    error('The input movie has not been background subtracted! Please perform background subtraction prior to spatial transformation!')    
end

%Check that all channels have been background subtracted
hasBS = movieData.processes_{iBSProc}.checkChannelOutput;   
if ~all(hasBS(p.ChannelIndex))
    error('Every channel selected for transformation must have been background subtracted! Please perform background subtraction first, or check the ChannelIndex parameter!')
end

%Check if bleedthrough correction has been performed on any channels.
iBTCProc = find(cellfun(@(x)(isa(x,'BleedthroughCorrectionProcess')),movieData.processes_),1);                          
%If so, check which channel(s).
if ~isempty(iBTCProc)
   hasBTC = movieData.processes_{iBTCProc}.checkChannelOutput;   
   if any(hasBTC)
      disp('Using bleed-through corrected images for channels:')
      arrayfun(@(x)(disp(num2str(x))),find(hasBTC));
   end
else
    hasBTC = false(1,numel(movieData.channels_));
end

nChanCorr = length(p.ChannelIndex);

%Set up the input /output directories for each channel
for j = 1:nChanCorr
    
    if hasBTC(p.ChannelIndex(j))
        %If available, use the bleed-through corrected images
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iBTCProc}.outFilePaths_{1,p.ChannelIndex(j)});
    else
        %Otherwise, use background subtracted
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iBSProc}.outFilePaths_{1,p.ChannelIndex(j)});
    end
    
    %The output is a sub-dir of the directory specified by OutputDirectory
    currDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];
    
    %Check/set up directory
    mkClrDir(currDir);
    
    movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex(j),currDir);                    
end

%Check if transform files have been specified, and if not, get them
if ~iscell(p.TransformFilePaths)
    %If only a single path was entered, it doesn't have to be a cella rray
    p.TransformFilePaths = {p.TransformFilePaths};
end

for j = 1:nChanCorr
    
    if isempty(p.TransformFilePaths{p.ChannelIndex(j)})
        [currFile,currPath] = uigetfile('*.mat',...
            ['Please select the transformation file for channel ' ...
            num2str(p.ChannelIndex(j)) ':']);
        
        if currFile == 0
            error('You must specify a transformation file to continue!')
        end        
        p.TransformFilePaths{p.ChannelIndex(j)} = [currPath currFile];                    
    end
    %This method will check validity of file....
    movieData.processes_{iProc}.setTransformFilePath(p.ChannelIndex(j),...
                                                p.TransformFilePaths{p.ChannelIndex(j)});        
end


%% ------- Init ------ %%



disp('Loading transformation...')
    


%Get the actual transformations for each channel
xForms = movieData.processes_{iProc}.getTransformation(p.ChannelIndex);
inNames = movieData.processes_{iProc}.getInImageFileNames(p.ChannelIndex);

%Get original image size. Image pixels that are transformed out of this
%area will be omitted to preserve this size
n = movieData.imSize_(1);
m = movieData.imSize_(2);        


%% ------- Spatial Transformation ------ %%
%Transform all images in requested channels and write them to a new
%directory.


disp('Transforming images....')

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, transforming correcting channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        

nImages = movieData.nFrames_;
nImTot = nImages * nChanCorr;

for iChan = 1:nChanCorr
    
    %Get directories for readability
    inDir  = movieData.processes_{iProc}.inFilePaths_{1,p.ChannelIndex(iChan)};    
    outDir = movieData.processes_{iProc}.outFilePaths_{1,p.ChannelIndex(iChan)};    
    
    disp(['Transforming images for channel ' num2str(p.ChannelIndex(iChan))])
    disp(['Transforming images from ' inDir ', results will be stored in ' outDir]);     
    disp(['Using transform file : ' p.TransformFilePaths{p.ChannelIndex(iChan)}]);
    
    if ~p.BatchMode        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, transforming channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end        



    for iImage = 1:nImages        
        
        currIm = imread([inDir filesep inNames{iChan}{iImage}]);                                
        
        currIm = imtransform(currIm,xForms{iChan},'XData',[1 m],'YData',[1 n],'FillValues',0);
        
        imwrite(currIm,[outDir filesep pString inNames{iChan}{iImage}]);               
        
        if ~p.BatchMode && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end                        
        
    end    
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end



%% ------- Mask Transformation ----- %%

if p.TransformMasks    
    
    %Get the indices of any previous mask intersection process
    iMaskTransfProc = movieData.getProcessIndex('MaskTransformationProcess',1,0);
    
    %If the process doesn't exist, create it
    if isempty(iMaskTransfProc)
        iMaskTransfProc = numel(movieData.processes_)+1;
        movieData.addProcess(MaskTransformationProcess(movieData,p.OutputDirectory));
    end
    maskTransfProc = movieData.processes_{iMaskTransfProc};
       
    %Set up the parameters for mask transformation
    maskTransfParams.ChannelIndex = p.ChannelIndex;
    maskTransfParams.TransformFilePaths = p.TransformFilePaths;    
    maskTransfParams.SegProcessIndex = p.SegProcessIndex;
    maskTransfParams.BatchMode = p.BatchMode;
    
    
    parseProcessParams(maskTransfProc,maskTransfParams);
    maskTransfProc.run;
        
end





%% ------ Output and Finalization ----- %%


%Store parameters/settings in movieData structure

movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished!')

