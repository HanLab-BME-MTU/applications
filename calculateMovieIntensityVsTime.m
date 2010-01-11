function movieData = calculateMovieIntensityVsTime(movieData,varargin)

% movieData = calculateMovieIntensityVsTime(movieData)
% 
% movieData = calculateMovieIntensityVsTime(movieData,'OptionName',optionValue...)
%
% This function calculates various statistics about how the intensity in
% the specified movie channels varies over time. 
% 
% Calculated values include:
% 
%   Average intensity, total intensity. Average masked intensity, total
%   masked intensity (if masks present).
% 
%
% Input:
%
%   movieData - The structure describing the movie to be analyzed, as
%               created by setupMovieData.m
%
%   Possible Option Names:
%       ('OptionName' -> possible values)
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) to perform intensity calculations on.
%       This index corresponds to the channel directories location in the
%       cell array movieData.channelDirectory. If not input, statistics
%       will be calculated for ALL channels
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
% 
% Output:
% 
%   The results will be saved to the movie's analysis directory and logged
%   in the movieData structure. Figures showing the intensity vs time
%   statistics will also be created and saved to the analysis folder.
%
% Hunter Elliott, 11/2009
% 
%% ------ Parameters ----- %%

pFileName = 'intensity_vs_time_'; % The string to prepend before the channel name when saving the results to file.


%% ------ Input ----- %%

movieData = setupMovieData(movieData,'intensityVsTime');

[batchMode,iChannels] = parseInput(varargin);

% --- Defaults --- %

if isempty(batchMode)
    batchMode = false;
end



if isempty(iChannels) %Default is to analyze all channels
    nChan = length(movieData.channelDirectory);
    iChannels = 1:nChan;
else
    nChan = length(iChannels);   
end


%% --------------- Init ---------- %%


allMeanI = cell(1,nChan);
allMeanMaskedI = cell(1,nChan);
allMeanBackgroundI = cell(1,nChan);

allTotalI = cell(1,nChan);
allTotalMaskedI = cell(1,nChan);
allTotalBackgroundI = cell(1,nChan);


hasMasks = false(1,nChan);
hasBackMasks = false(1,nChan);



%% ----- Calculate Intensity vs time in all Channels -----%%

disp('Starting intensity calculations...')

for iChan = 1:nChan

    nImages = movieData.nImages(iChannels(iChan));
    allMeanI{iChan} = zeros(1,nImages);
    allTotalI{iChan} = zeros(1,nImages);
    
    imageDir = [movieData.imageDirectory filesep movieData.channelDirectory{iChannels(iChan)}];
    imageFileNames = dir([imageDir filesep '*.tif']);    
            
    
    
    hasMasks(iChan) = checkMovieMasks(movieData,iChannels(iChan));   
    if hasMasks(iChan)
        maskDir = [movieData.masks.directory filesep movieData.masks.channelDirectory{iChannels(iChan)}];
        maskFileNames = dir([maskDir filesep '*.tif']);
        allMeanMaskedI{iChan} = zeros(1,nImages);
        allTotalMaskedI{iChan} = zeros(1,nImages);
    end
    
    %Check background masks for this channel
    hasBackMasks(iChan) = checkMovieBackgroundMasks(movieData,iChannels(iChan));    
    if hasBackMasks(iChan)
        backMaskDir = [movieData.backgroundMasks.directory ...
            filesep movieData.backgroundMasks.channelDirectory{iChannels(iChan)}];
        backMaskFileNames = dir([backMaskDir filesep '*.tif']);
        allMeanBackgroundI{iChan} = zeros(1,nImages);
        allTotalBackgroundI{iChan} = zeros(1,nImages);
    end
    
    disp(['Calculating intensity vs time for channel ' movieData.channelDirectory{iChannels(iChan)}])
    
    for iImage = 1:nImages
        
        currIm = imread([imageDir filesep imageFileNames(iImage).name]);
        
        if hasMasks(iChan)
            currMask = imread([maskDir filesep maskFileNames(iImage).name]);
            allMeanMaskedI{iChan}(iImage) = mean(currIm(currMask(:)));    
            allTotalMaskedI{iChan}(iImage) = sum(currIm(currMask(:)));    
        end
        
        if hasBackMasks(iChan)
            currMask = imread([backMaskDir filesep backMaskFileNames(iImage).name]);
            allMeanBackgroundI{iChan}(iImage) = mean(currIm(currMask(:)));
            allTotalBackgroundI{iChan}(iImage) = sum(currIm(currMask(:)));
        end
        
        allMeanI{iChan}(iImage) = mean(currIm(:));                    
        allTotalI{iChan}(iImage) = sum(currIm(:));
        
    end
    
end


%% ------- Output -------%%



%----- Make and save Figures -----%


% Whole-Image Stats Figure

disp('Making figures...')

if batchMode
    figHan = figure('Visible','off');
else
    figHan = figure;
end

chanColors = jet;
chanColors = chanColors(round(linspace(1,length(chanColors),nChan)),:);

subplot(2,1,1)
hold on
tData = 0:max(movieData.nImages) * movieData.timeInterval_s;
for i = 1:nChan
    plot(tData(1:length(allMeanI{i})),allMeanI{i},'color',chanColors(i,:))
end
set(gca,'Color','None')
ylabel('Average Intensity')
legend(movieData.channelDirectory(iChannels),'Interpreter','none')

subplot(2,1,2)
hold on
for i = 1:nChan
    plot(tData(1:length(allMeanI{i})),allTotalI{i},'color',chanColors(i,:))
end
ylabel('Total Intensity')
xlabel('Time, seconds')
set(gca,'Color','None')


hgsave(figHan,[movieData.intensityVsTime.directory filesep pFileName '.fig'])



% Masked Stats Figure

if any(hasMasks)
    if batchMode
        figHanM = figure('Visible','off');
    else
        figHanM = figure;
    end
    subplot(2,1,1)
    hold on    
    for i = 1:nChan
        if hasMasks(i)
            plot(tData(1:length(allMeanMaskedI{i})),allMeanMaskedI{i},'color',chanColors(i,:))
        end
    end
    set(gca,'Color','None')
    ylabel('Average Masked Intensity')
    legend(movieData.masks.channelDirectory(iChannels(hasMasks)),'Interpreter','none')

    subplot(2,1,2)
    hold on
    for i = 1:nChan
        if hasMasks(i)
            plot(tData(1:length(allTotalMaskedI{i})),allTotalMaskedI{i},'color',chanColors(i,:))
        end
    end
    ylabel('Total Masked Intensity')
    xlabel('Time, seconds')
    set(gca,'Color','None')
    hgsave(figHanM,[movieData.intensityVsTime.directory filesep 'masked_' pFileName '.fig'])

end

% Background Stats Figure

if any(hasBackMasks)
    if batchMode
        figHanBM = figure('Visible','off');
    else
        figHanBM = figure;
    end
    subplot(2,1,1)
    hold on    
    for i = 1:nChan
        if hasBackMasks(i)
            plot(tData(1:length(allMeanBackgroundI{i})),allMeanBackgroundI{i},'color',chanColors(i,:))
        end
    end
    set(gca,'Color','None')
    ylabel('Average Background Intensity')
    legend(movieData.masks.channelDirectory(iChannels(hasBackMasks)),'Interpreter','none')

    subplot(2,1,2)
    hold on
    for i = 1:nChan
        if hasBackMasks(i)
            plot(tData(1:length(allTotalBackgroundI{i})),allTotalBackgroundI{i},'color',chanColors(i,:))
        end
    end
    ylabel('Total Background Intensity')
    xlabel('Time, seconds')
    set(gca,'Color','None')

    hgsave(figHanBM,[movieData.intensityVsTime.directory filesep 'background_' pFileName '.fig'])

end
    

% ----- write intensity stats to file ---- %


for i = 1:nChan %Each channel is written to a seperate file
       
    meanIntensity = allMeanI{i}; %Seperate the channels for saving
    totalIntensity = allTotalI{i};
    
    saveString = {'meanIntensity','totalIntensity'};
    
    if hasBackMasks(i)           
        meanBackgroundIntensity = allMeanBackgroundI{i};
        totalBackgroundIntensity = allTotalBackgroundI{i};
        saveString = [saveString {'meanBackgroundIntensity' 'totalBackgroundIntensity'}];                        
    end
    if hasMasks(i)
        meanMaskedIntensity = allMeanMaskedI{i};
        totalMaskedIntensity = allTotalMaskedI{i};
        saveString = [saveString {'meanMaskedIntensity','totalMaskedIntensity'}];
    end
    
    save([movieData.intensityVsTime.directory filesep ...
        pFileName movieData.channelDirectory{iChannels(i)} '.mat'],saveString{:})
        
end

movieData.intensityVsTime.status = 1;
movieData.intensityVsTime.channels(iChannels) = true;
movieData.intensityVsTime.maskedChannels(iChannels) = hasMasks;
movieData.intensityVsTime.backgroundMaskedChannels(iChannels) = hasBackMasks;
movieData.intensityVsTime.dateTime = datestr(now);
movieData.intensityVsTime.fileNamePrefix = pFileName;

updateMovieData(movieData);

disp('Finished!')

function [batchMode, iChannels] = parseInput(argArray)

%Init output
batchMode = [];
iChannels = [];

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
    
   switch argArray{i}                     
              
       case 'BatchMode'
           batchMode = argArray{i+1};
           
       case 'ChannelIndex'
           iChannels = argArray{i+1};
                      
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end


