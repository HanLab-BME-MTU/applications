function movieData = sampleMovieWindows(movieData,paramsIn)
%SAMPLEMOVIEWINDOWS samples the movie images in the sampling windows created with getMovieWindows.m 
% 
% movieData = sampleMovieWindows(movieData)
% movieData = sampleMovieWindows(movieData,paramsIn)
% 
% This function goes through each frame in the movie and samples the images
% in the areas occupied by each sampling window that has been created using
% getMovieWindows.m. Various statistics for the pixels inside each window
% are calculated, and stored in an array the same size as the window array. 
% 
% 
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       samples to. Samples for different channels will be saved as
%       files in this directory.
%       If not input, the samples will be saved to the same directory as the
%       movieData, in a sub-directory called "window_samples"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channels to sample images from. If not input,
%       all available channels will be sampled.
%
%       ('ProcessIndex' -> Positive integer scalar) Optional. This
%       specifies the index of an ImageProcessingProcess or
%       DoubleProcessingProcess in the movieData's process array to use the
%       output images of for sampling. If not specified, the raw images
%       will be sampled.
%
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - The updated MovieData object, with the parameters and
%   locations of the samples stored in it.
% 
%   Additionally, the window samples for each channel will be written to a
%   file in the OutputDirectory, as .mat files, with one file per channel.
% 
% 
% Hunter Elliott
% 7/2010
%
%% --------- Parameters ------------ %%

pString = 'window_samples_for_channel_'; %Prefix for saving samples to file

%% ---------- Input ---------------- %%


if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];  
end
   
%Check if the movie has been sampled before
iProc = movieData.getProcessIndex('WindowSamplingProcess',1,false);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(WindowSamplingProcess(movieData,movieData.outputDirectory_));
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has been windowed, and find the desired process.
iWinProc = movieData.getProcessIndex('WindowingProcess',1,~p.BatchMode);
if isempty(iWinProc)
    error('The movie could not be sampled, because it has not been windowed yet!')
end

%Make sure that the windows are okay.
if ~movieData.processes_{iWinProc}.checkChannelOutput;
    error('The window files for the input movie are not valid!')    
end



%% -------- Init ---------- %%


nChan = numel(p.ChannelIndex);
nFrames = movieData.nFrames_;
imSize = movieData.imSize_;


%Get & set the input image directories and file names
if isempty(p.ProcessIndex)
    
    imNames = movieData.getImageFileNames(p.ChannelIndex);
    imDirs = movieData.getChannelPaths(p.ChannelIndex);
    isDouble = false;
    
else
   
    %Make sure the process specified is an ImageProcessingProcess
    if isa(movieData.processes_{p.ProcessIndex},'DoubleProcessingProcess')
        isDouble = true;
    elseif isa(movieData.processes_{p.ProcessIndex},'ImageProcessingProcess')
        isDouble = false;        
    else
        error('The process selected for input by the ProcessIndex parameter must be an ImageProcessingProcess or a DoubleProcessingProcess!');
    end
    
    imNames = movieData.processes_{p.ProcessIndex}.getOutImageFileNames(p.ChannelIndex);
    imDirs = movieData.processes_{p.ProcessIndex}.outFilePaths_(p.ChannelIndex);

end


%Set up and store the output directories for the window samples.
mkClrDir(p.OutputDirectory)
    
%Initialize sample array
samples(1:nChan) = struct('avg',[],'std',[],'max',[],'min',[],'med',[]);
fNames = fieldnames(samples(1));
nFields = numel(fNames);

for j = 1:nFields

    for k = 1:nChan
        
        samples(k).(fNames{j}) = nan(...
                  movieData.processes_{iWinProc}.nSliceMax_,...
                  movieData.processes_{iWinProc}.nBandMax_,...
                  nFrames);
    end
end


%Get the mask information from the windowing process
iSegProc = movieData.processes_{iWinProc}.funParams_.SegProcessIndex;
if ~isa(movieData.processes_{iSegProc},'SegmentationProcess')
    error('The segmentation process specified by the windowing process is invalid! Please check settings and re-run windowing!')
end

%Store these in the parameter structure.
p.SegProcessIndex = iSegProc;
p.MaskChannelIndex = movieData.processes_{iWinProc}.funParams_.ChannelIndex;
nMaskChan = numel(p.MaskChannelIndex);

%Get the mask directories and file names
maskDir = movieData.processes_{iSegProc}.outFilePaths_(p.MaskChannelIndex);
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.MaskChannelIndex);



%% --------- Sampling --------- %%

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, sampling windows...');
end  

disp('Starting window sampling...');


for iFrame = 1:nFrames
     
    %Load the windows
    currWin = movieData.processes_{iWinProc}.loadChannelOutput(iFrame);    
        
    %Load the mask(s) to use first, so we can combine them and use this to
    %mask every channel.
    currMask = true(imSize);
    for j = 1:nMaskChan
        currMask = currMask & imread([maskDir{j} filesep maskNames{j}{iFrame}]);
    end    
    
    %Go through each channel and sample it
    for iChan = 1:nChan
               
        if ~isDouble
            currIm = imread([imDirs{iChan} filesep imNames{iChan}{iFrame}]);                                        
        else
            currIm = movieData.processes_{p.ProcessIndex}.loadOutImage(p.ChannelIndex(iChan),iFrame);
        end
        
        currSamples = sampleImageWindows(currWin,currIm,currMask);
        
        %Copy these into the whole-movie array
        currSize = size(currSamples.(fNames{1}));%all field arrays are same size
        for j = 1:nFields            
            samples(iChan).(fNames{j})(1:currSize(1),1:currSize(2),iFrame) ...
                = currSamples.(fNames{j});                         
        end
        
        
    end                    
        
    if ~p.BatchMode && mod(iFrame,5)
        %Update the waitbar occasionally  to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
    
end

%% ------- Output ------- %%

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar);
end

for j = 1:nChan
    %Save the samples to file    
    actSamples = samples(j); %#ok<NASGU>
    fPath = [p.OutputDirectory filesep pString num2str(p.ChannelIndex(j)) '.mat'];
    save(fPath,'actSamples');
    movieData.processes_{iProc}.setOutFilePath(p.ChannelIndex(j),fPath);    
end

%Update the movie data, save it
movieData.processes_{iProc}.setDateTime;
movieData.processes_{iProc}.setPara(p);%We've stored additional parameters, so add to the process structure.
movieData.save; %Save the new movieData to disk


disp('Finished sampling!')




