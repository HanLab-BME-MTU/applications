function movieData = outputMovieRatios(movieData,paramsIn)
%OUTPUTMOVIERATIOS writes the ratio images to .tif files
% 
% movieData = outputMovieRatios(movieData)
% 
% movieData = outputMovieRatios(movieData,paramsIn)
%
% This function converts the ratio images for a selected channel (which are
% floating point, double-precision images) into integer-valued .tif images,
% and saves them to a specified output directory.
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
%       string specifying the directory to output the ratio images to. If
%       not input, the corrected images will be saved to the same directory
%       as the movieData, in a sub-directory called "ratio_tiffs"
%
%       ('ChannelIndex'-> Positive integer scalar) The integer index of the
%       NUMERATOR of the ratio channel to perform photbleach correction on.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, the user will be asked to select
%       from the movie's channels.
%
%       ('ScaleFactor'-> Positive scalar)
%       This specifies the value to multiply the original ratios by when
%       converting them to .tif images. It is recommended that this number
%       be at least 1000 or higher to avoid rounding error in the resulting
%       .tif images.
%
%       ('BatchMode' -> True/False) If true, all graphical outputs and user
%       interaction is suppressed. 
% 
%
%
% Output:
%
%   movieData - the updated movieData object with the correction
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The ouput images are written to the directory specified by the
%   parameter OuptuDirectory, in a sub-directory. They will be stored as
%   bit-packed .tif files. 
%
% 
% Hunter Elliott
% 6/2010
%

%% ------ Input ------ %%

%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Make sure the movie has been ratioed
iRProc = find(cellfun(@(x)(isa(x,'RatioProcess')),movieData.processes_),1);                          

if isempty(iRProc)
    error('The input movie has not been ratioed! Please perform ratioing prior to outputing ratio images!')
end

%Look for previous output processes
iProc = find(cellfun(@(x)(isa(x,'OutputRatioProcess')),movieData.processes_),1);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(OutputRatioProcess(movieData,movieData.outputDirectory_));                                                                                                 
end


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


if isempty(p.ChannelIndex)
    if ~p.BatchMode
        p.ChannelIndex = selectMovieChannels(movieData,1,'Select the ratio channel to output:');
    else
        error('In batch mode, you must specify the channel to output!')
    end
end

%Check if the movie has a photbleach correction processes 
iPBProc = find(cellfun(@(x)(isa(x,'PhotobleachCorrectionProcess')),movieData.processes_),1);                          

if isempty(iPBProc)
    hasPB = false;
else
    hasPB = movieData.processes_{iPBProc}.checkChannelOutput(p.ChannelIndex);
end

if length(p.ChannelIndex) ~=1
    error('You can only output one ratio channel at a time!')
end
%% -------- Init -------- %%

disp('Starting ratio output ...')


%Set up input and output directories
if hasPB
    inDir = movieData.processes_{iPBProc}.outImagePaths_{p.ChannelIndex};
    inNames = movieData.processes_{iPBProc}.getOutImageFileNames(p.ChannelIndex);
else
    inDir = movieData.processes_{iRProc}.outImagePaths_{p.ChannelIndex};
    inNames = movieData.processes_{iRProc}.getOutImageFileNames(p.ChannelIndex);        
end

outDir = p.OutputDirectory;
mkClrDir(outDir);
movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex,outDir);


nImages = movieData.nFrames_;

%% ------- Output -------- %%
%Writes the ratio images to file.

disp(['Outputing ratio images for channel ' inDir ])
disp(['Writing ratio .tiff images to folder ' outDir])

%Disable convert-to-integer warning
warning('off','MATLAB:intConvertNonIntVal');

if ~p.BatchMode        
    wtBar = waitbar(0,'Please wait, writing ratios to .tif images...');
end        


%Go through all the images and correct them
for iImage = 1:nImages
   
    %Load the image
    currRat = load([inDir filesep inNames{1}{iImage}]);   
    fNames = fieldnames(currRat);
    currRat = currRat.(fNames{1});
    
    %Scale the image
    currRat = currRat .* p.ScaleFactor;
    
    %Write it back to file.    
    imwrite(uint16(currRat),[outDir filesep inNames{1}{iImage}(1:end-4) '.tif'],...
        'Compression','none'); %Disable compression to ensure compatibility
    
    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iImage / nImages,wtBar)
    end                        

end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

%% ----- Finalization ---- %%


%Log the output in the movieData object and save it

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData;


disp('Finished!')




