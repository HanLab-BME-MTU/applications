function [tracksNA] = getTrackSignals(MD,tracksNA,iChan,preFrames)
% [] = getTrackSignals(MD,tracksNA,iChan)
% obtain the signal in iChan channel of movie in movieData MD. 
% 
% input:    MD:    movieData file 
%           outputPath                  outputPath
%           tracksNA            tracks
%           iChan               channel number of interest
%           preFrames       the number of frames that will be used to
%           collect the signal ahead of starting time point of each track
% output:  tracksNA          tracks that contains intensity information in
%                                           ampChan

% Sangyoon Han July 2014

if nargin <4
    preFrames = 0;
end
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
    iStart = max(1, tracksNA(k).startingFrame-preFrames);
    iEnd = tracksNA(k).endingFrame;
    for ii=iStart:iEnd
        % Get the intensity value from tracks
        if ~tracksNA(k).presence(ii) && iStart >= tracksNA(k).startingFrame
            continue
        end
        if ii<tracksNA(k).startingFrame
            p = tracksNA(k).startingFrame;
        else
            p = ii;
        end
        curX = tracksNA(k).xCoord(p);
        curY = tracksNA(k).yCoord(p);
        tracksNA(k).ampChan(ii) = ichanImg(round(curY),round(curX),ii);
    end
end


