function movieData = processBleedthroughMovie(movieData,fluorChannel,bleedChannel)
%PROCESSBLEEDTHROUGHMOVIE processes and extracts bleedthrough coeefficients of the input movie
% 
% movieData = processBleedthroughMovie(movieData)
% 
% movieData = processBleedthroughMovie(movieData,fluorChannel,bleedChannel)
% 
% This function is for processing bleedthrough correction movies - that is,
% movies where a single fluorophore was imaged in two channels to determine
% how much that fluorophore "bleeds through" into the channel used to image
% another fluorophore. The movie is first processed by performing
% segmentation, dark-current subtraction, background subtraction, and
% shading correction and then the bleedthrough coefficients are calculated
% using calculateMovieBleedthrough.m
% 
%
% Input:
% 
%
%   movieData - the structure describing the movie, as created with
%   setupMovieData.m
%
%   fluorChannel - A character array containing the name of the channel
%   which has a fluorophore.
%
%   bleedChannel - A character array containing the name of the channel
%   which does NOT have a fluorohpore - the "bleed through" channel.
% 
% 
% Output:
%
%
%   movieData - the modified structure with the bleedthrough calculation
%   logged in it. 
%
%
% Hunter Elliott 
% 2/2010
%

%% ----------- Input --------- %%

if nargin < 1
    movieData = [];
end

movieData = setupMovieData(movieData);

if nargin < 2 || isempty(fluorChannel);
    fluorChannel = selectMovieChannels(movieData,0,'Please select the channel WITH a fluorophore:');
else
    %Get the index from the name
    fluorChannel = find(strcmp(fluorChannel,movieData.channelDirectory),1);    
end
if nargin < 3 || isempty(bleedChannel)
    bleedChannel = selectMovieChannels(movieData,0,'Please select the channel WITHOUT a fluorophore:');
else
    bleedChannel = find(strcmp(bleedChannel,movieData.channelDirectory),1);
end

%% -------- Init ------- %%

%create a regexp to search for correctly named
%channels for the shade correction images.
nameFshadeCorr = ['(' movieData.channelDirectory{fluorChannel} ...
    '.*shade)|(shade.*' movieData.channelDirectory{fluorChannel} ')'];
nameBshadeCorr = ['(' movieData.channelDirectory{bleedChannel} ...
    '.*shade)|(shade.*' movieData.channelDirectory{bleedChannel} ')'];

%Look for appropriately named channels for the shade correction images
iSCfluor = find(cellfun(@(x)(~isempty(regexpi(x,nameFshadeCorr))),movieData.channelDirectory));
iSCbleed = find(cellfun(@(x)(~isempty(regexpi(x,nameBshadeCorr))),movieData.channelDirectory));

if length(iSCfluor) ~= 1 %If it failed, or was ambiguous, ask the user
   iSCfluor = selectMovieChannels(movieData,0,...
       ['Please select the shade correction channel for channel ' movieData.channelDirectory{fluorChannel}]);
end
if length(iSCbleed) ~= 1 %If it failed, or was ambiguous, ask the user
   iSCbleed = selectMovieChannels(movieData,0,...
       ['Please select the shade correction channel for channel ' movieData.channelDirectory{fluorChannel}]);
end


%% ---------- Processing ------ %%

%Create Masks for the movie - we only segment the fluorophore channel - it
%has much better SNR, and it's masks will be used for both channels.
movieData = thresholdMovie(movieData,'ChannelIndex',fluorChannel);

%Create background masks
movieData = createMovieBackgroundMasks(movieData,'ChannelIndex',fluorChannel);

%Apply dark-current correction
movieData = darkCurrentCorrectMovie(movieData,'ChannelIndex',[fluorChannel bleedChannel]);
newFchan = find(movieData.darkCurrentCorrection.iFrom == fluorChannel); %Find the indices for the dark-current corrected images
newBchan = find(movieData.darkCurrentCorrection.iFrom == bleedChannel);

%Apply shade correction
movieData = shadeCorrectMovie(movieData,'ChannelIndex',[newFchan newBchan],...
    'ShadeImageChannels',[iSCfluor iSCbleed]);
    
%Apply background subtraction to the shade-corrected images
newFchan = find(movieData.shadeCorrection.iFrom == newFchan); %Find the indices for the shade-corrected images
newBchan = find(movieData.shadeCorrection.iFrom == newBchan);
movieData = backgroundSubtractMovie(movieData,'ChannelIndex',[newFchan newBchan],...
    'BackgroundMaskIndex',[fluorChannel fluorChannel]);%and use the masks from the orginal fluorophore channel, since it has much better SNR

%Get indices of the background-subtracted images
newFchan = find(movieData.backgroundSubtraction.iFrom == newFchan);
newBchan = find(movieData.backgroundSubtraction.iFrom == newBchan);

%Calculate bleedthrough coefficients using these processed images
movieData = calculateMovieBleedthrough(movieData,...
    'FluorophoreChannel',newFchan,'BleedthroughChannel',newBchan,...
    'FluorophoreMaskChannel',fluorChannel,'BleedthroughMaskChannel',fluorChannel);

    
    
    