function [tracksNA] = getTrackSignals(MD,tracksNA,iChan)
% [] = getTrackSignals(MD,tracksNA,iChan)
% obtain the signal in iChan channel of movie in movieData MD. 

% input:    MD:    movieData file 
%           outputPath                  outputPath
%           tracksNA            tracks
%           iChan               channel number of interest
% output:  tracksNA          tracks that contains intensity information in
%                                           ampChan

% Sangyoon Han July 2014

%% Data Set up
% Get whole frame number
nFrames = MD.nFrames_;
% % name of the intensity of iChan
% nameChan = ['amp' num2str(iChan)];
%% Get the image stack
ichanImg = zeros(MD.imSize_(1),MD.imSize_(2),nFrames);
for ii=1:nFrames
    % Get the image in iChan
    ichanImg(:,:,ii) =  MD.getChannel(iChan).loadImage(ii);
end
%% Get the image in iChan
for k=1:numel(tracksNA)
    iStart = tracksNA(k).startingFrame;
    iEnd = tracksNA(k).endingFrame;
    for ii=iStart:iEnd
        % Get the intensity value from tracks
        if ~tracksNA(k).presence(ii)
            continue
        end
        curX = tracksNA(k).xCoord(ii);
        curY = tracksNA(k).yCoord(ii);
        tracksNA(k).ampChan(ii) = ichanImg(round(curY),round(curX),ii);
    end
end


