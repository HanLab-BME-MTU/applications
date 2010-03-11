function movieData = processBiosensorMovie(movieData,varargin)
%PROCESSBIOSENSORMOVIE produces ratio maps from raw biosensor imaging data 
% 
% movieData = processBiosensorMovie(movieData)
% 
% movieData = processBiosensorMovie(movieData,'OptionName',optionValue)
%
% This function produces ratio images from the input biosensor movie by
% performing a series of processing steps. If some processing steps have
% been performed, only those which have yet to be completed will be
% performed. 
% 
%   The processing steps, and their order, are as follows:
%
%       1 - Mask creation
%       2 - Background mask creation
%       3 - Mask Refinement (Optional - not run by default)
%       4 - Dark-Current correction 
%       5 - Shade correction
%       6 - Background Subtraction
%       7 - Bleedthrough correction (Optional - not run by default)
%       8 - Transformation (Optional-not run by default)
%       9 - Ratioing
%      10 - Photobleach correction
% 
% Input:
% 
%   movieData - movie information structure, as created by setupMovieData.m
%   Optional. If not input, the user will be asked to create a movieData by
%   specifying an analysis directory and image directory.
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%
%       ('VolumeChannel' -> Positive Integer) The integer index of the
%       volume-correction channel. (The denominator of the ratio). This
%       number corresponds to this channel's index in the field
%       movieData.channelDirectory.
%
%       ('ActivityChannel' -> Positive Integer) The integer index of the
%       activity channel. (The numerator of the ratio). This number
%       corresponds to this channel's index in the field
%       movieData.channelDirectory.
%
%       ('RunSteps' -> integer vector) This vector indicates which steps to
%       run. It should be nSteps long. 1 Indicates the step should be run
%       no matter what. 0 Indicates the step should only be run if it
%       hasn't been run yet
%      -1 Indicates the step should NOT be run no matter what.
%       For example, [0 0 -1 1 0 0 0 0] Means that step 3 should be
%       skipped, step four should be run no matter what, and all other
%       steps should only be run if they haven't been completed previously.
%       Optional. Default is to run all uncompleted steps, except the
%       transformation & mask refinemtn which are skipped by default.
%
%       ('BatchMode' -> True/False) If this option value is set to true,
%       all graphical output and user interaction is suppressed. (No
%       progress bars, no dialogue boxes)
%       Optional. Default is False.
%
%   Processing Step Options:
%   
%       Options for individual steps can be passed to the functions
%       performing each step as follows:
%
%       ('processName_OptionName' -> optionValue)
%       Where processName is the name of the processing step, such as
%       "masks" or "transformation", OptionName is the name of the option
%       in the processing function, and optionValue is the value to pass.
%       For example:
%
%       processBiosensorMovie(movieData,'backGroundMasks_ChannelIndex', [1
%       3]) Would cause the background mask creation step to use channels
%       number 1 and 3.
%
% Output: 
% 
%   The resulting ratio images, and the images resulting from each
%   processing step will be written to new channel directories in the
%   movie's image directory. All processing steps and all parameters used
%   will be logged in the movieData structure.
%
% Hunter Elliott
% 11/2009
%

%% ------ Parameters ------ %%

%Array of function handles to use at each step in the processing.
stepFunctions = {...
    @thresholdMovie,...    
    @createMovieBackgroundMasks,...
    @refineMovieMasks,...
    @darkCurrentCorrectMovie,...
    @shadeCorrectMovie,...
    @backgroundSubtractMovie,...
    @bleedthroughCorrectMovie,...
    @transformMovie,...
    @ratioMovie,...        
    @photobleachCorrectMovieRatios,...
    };
    
%The name of the procedure at each step. Also the same as the name of the
%field in the movieData structure where each step is logged.
%(This is seperate because two functions may perform the same procedure)
procName = {...
    'masks',...
    'backgroundMasks',... 
    'maskRefinement',... 
    'darkCurrentCorrection',...
    'shadeCorrection',...
    'backgroundSubtraction',...
    'bleedthroughCorrection',...
    'transformation',...
    'ratioing',...    
    'photobleachCorrection',...
    };



nSteps = length(stepFunctions); %Number of steps in processing.

if nSteps ~=length(procName)%Make sure someone didn't fuck something up.
    error('Number of processing steps and function handles not equivalent! Be careful when editing this file!!!!')
end

%% ------ Input -------%%


if nargin < 1
    movieData = [];
end

movieData = setupMovieData(movieData,'biosensorProcessing');

[batchMode,runSteps,stepOptions,iAct,iVol] = parseInput(varargin,procName);
    
%---Defaults----%

if isempty(batchMode)
    batchMode = false;
end
if isempty(runSteps)
    runSteps = zeros(nSteps,1); %Default is to run only uncompleted steps...
    runSteps(strcmp('transformation',procName)) = -1;%Except the transformation...
    runSteps(strcmp('bleedthroughCorrection',procName)) = -1;%...bleedthrough correction...
    runSteps(strcmp('maskRefinement',procName)) = -1;%... and mask refinement, which are skipped by default.
end

if isempty(iVol)
    if ~batchMode
        iVol = selectMovieChannels(movieData,0,'Please select the volume correction channel:');
    else
        error('In batch mode you must specify a volume channel index!')
    end
end
if isempty(iAct)
    if ~batchMode
        iAct = selectMovieChannels(movieData,0,'Please select the activity channel:');
    else
        error('In batch mode you must specify a volume channel index!')
    end
end

nChanTot = length(movieData.channelDirectory);
if iAct > nChanTot || iAct < 1 || round(iAct) ~= iAct
    error('Invalid activity channel index! Check activity channel!')
end
if iVol > nChanTot || iVol < 1 || round(iVol) ~= iVol
    error('Invalid volume channel index! Check activity channel!')
end



%% ------ Init ------ %%

%Always use the background masks from the original images
iBack = find(strcmp('backgroundSubtraction',procName));%Determine which step is background subtraction
if isempty(stepOptions{iBack})
    stepOptions{iBack} = {'BackgroundMaskIndex',[iAct iVol]};
elseif ~any(strcmp('BackgroundMaskIndex',stepOptions{iBack})) %Check if the user already specified masks
    stepOptions{iBack} = [stepOptions{iBack} {'BackgroundMaskIndex',[iAct iVol]}];%If not, specify them
end

%Determine which step is the ratioing, because after this there is one
%channel instead of two
iRatioStep = find(strcmp('ratioing',procName));
%Use foreground masks from the original images for the ratio
if isempty(stepOptions{iRatioStep})
    stepOptions{iRatioStep} = {'MaskChannelIndex', [iAct iVol]};
elseif ~any(strcmp('MaskChannelIndex',stepOptions{iBack}))
    stepOptions{iRatioStep} = [stepOptions{iRat} {'MaskChannelIndex', [iAct iVol]}];
end



%% ----- Processing ----- %%
%Go through each step and call the processing function.


if ~batchMode
    wtBar = waitbar(0,'Please wait, processing biosensor images...');
end

iInput = [iAct iVol];

for iStep = 1:nSteps                               

    if runSteps(iStep) == 1 || (runSteps(iStep) == 0 ...
            && ~checkMovieProcedure(movieData,procName{iStep},iInput))
        
        if ~isempty(stepOptions{iStep})        
            argIn = ['ChannelIndex',iInput,'BatchMode',batchMode,stepOptions{iStep}];
            movieData = stepFunctions{iStep}(movieData,argIn{:});
        else
            movieData = stepFunctions{iStep}(movieData,'ChannelIndex',iInput,...
            'BatchMode',batchMode);
        end    

    end
    if ~(runSteps(iStep) == -1) && ~(iStep == nSteps) 
        %Unless the step was skipped, or it is the last step, we need to
        %get the output from the last step        
        if iStep < iRatioStep

            iOldIn = iInput;
            iInput(1) = find(movieData.(procName{iStep}).iFrom == iOldIn(1));
            iInput(2) = find(movieData.(procName{iStep}).iFrom == iOldIn(2));

        elseif strcmp(procName{iStep},'ratioing')
            iInput = find(movieData.(procName{iStep}).iFrom(1,:) == iInput(1));           
        else
            iInput = find(movieData.(procName{iStep}).iFrom == iInput(1));           
        end
    end
    if ~batchMode
        waitbar(iStep/nSteps,wtBar)
    end
end


%% ----- Output ----- %%

movieData.biosensorProcessing.status = 1;
movieData.biosensorProcessing.dateTime = datestr(now);
movieData.biosensorProcessing.iVolume = iVol;
movieData.biosensorProcessing.iActivity = iAct;

updateMovieData(movieData);



if ~batchMode && ishandle(wtBar)
    close(wtBar);
end


function [batchMode,runSteps,stepOptions,iAct,iVol] = parseInput(argArray,procNames)

nSteps = length(procNames);

%Init output
batchMode = [];
runSteps = [];
stepOptions = cell(nSteps,1);
iVol = [];
iAct = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg

    iUnderscore = regexp(argArray{i},'_','ONCE');%Check if this is an option for a particular step    
    
    if ~isempty(iUnderscore)
        %Figure out which step this option is for
        sNum = strcmp(argArray{i}(1:iUnderscore-1),procNames);
        
        if isempty(stepOptions{sNum}) %Get the option name and value, removing the step number            
            stepOptions{sNum} = {argArray{i}(iUnderscore+1:end) argArray{i+1}}; 
        else            
            stepOptions{sNum} = [stepOptions{sNum} {argArray{i}(iUnderscore+1) argArray{i+1}}];                              
        end                
    else
        switch argArray{i}                     



           case 'BatchMode'
               batchMode = argArray{i+1};

           case 'RunSteps'
               runSteps = argArray{i+1};

           case 'VolumeChannel'

               iVol = argArray{i+1};

           case 'ActivityChannel'

               iAct = argArray{i+1};

           otherwise

               error(['"' argArray{i} '" is not a valid option name! Please check input!'])
        end
    end
end


