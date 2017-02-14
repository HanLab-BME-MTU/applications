function movieInfoLastFrame = genMovieInfoLastFrameFromTracksSparse(tracksSim)
%GENMOVIEINFOFROMTRACKS generates a list of detected features per frame from supplied tracks
%
%SYNOPSIS movieInfoLastFrame = genMovieInfoFromTracksTracksSparse(tracksSim,percentMissing)
%
%INPUT  tracksSim     : simulate tracks (compTracks).
%      
%OUTPUT movieInfoLastFrame: List of detected features in the last frame, in the format
%                  required for the input of trackWithGapClosing and
%                  trackCloseGapsKalman.
%
%                  .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                  .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             
%                  .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%
%
%Khuloud Jaqaman, October 2007
% Modified by Luciana de Oliveira, January 2017.

%call imput variables from structure
compTracks=tracksSim.compTracks;

%get number of frames in movie
seqOfEvents = vertcat(compTracks.seqOfEvents);
lastFrame = max(seqOfEvents(:,1));


%define the start time 
startTime = seqOfEvents(1,1);
%initialize variables
    xCoord= [];
    yCoord= [];
    amp   = [];

%transform each track in compTrack in the full matrix version
 
for     i=1:length(compTracks); 
    compTracksFull.tracksFeatIndxCG = full(compTracks(i).tracksFeatIndxCG);
    compTracksFull.tracksCoordAmpCG = full(compTracks(i).tracksCoordAmpCG); 
    compTracksFull.seqOfEvents = compTracks(i).seqOfEvents; 
    [trackedFeatureInfo,~,~,~,~] = convStruct2MatIgnoreMS(compTracksFull);
    
% take coordinates and intensity     
     xCoord = [xCoord; ...
                trackedFeatureInfo(:,(lastFrame-startTime)*8+1),trackedFeatureInfo(:,(lastFrame-startTime)*8+5)];
     yCoord = [yCoord; ...
                trackedFeatureInfo(:,(lastFrame-startTime)*8+2),trackedFeatureInfo(:,(lastFrame-startTime)*8+6)];
     amp    = [amp; ...
                trackedFeatureInfo(:,(lastFrame-startTime)*8+4),trackedFeatureInfo(:,(lastFrame-startTime)*8+8)];

end

% remove the tracks that are not present in the last frame (all row will be zero)
  
xCoord=xCoord(any(xCoord,2),:);
yCoord=yCoord(any(yCoord,2),:);
amp=amp(any(amp,2),:);
    


%store feature information in movieInfoLastFrame
    movieInfoLastFrame.xCoord = xCoord;
    movieInfoLastFrame.yCoord = yCoord;
    movieInfoLastFrame.amp = amp;
    

%%% ~~ the end ~~ %%%
