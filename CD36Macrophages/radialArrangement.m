function [radialityParam,centerCoord,numTracksLin] = radialArrangement(tracks,...
    centerCoordInit,minTrackLen,diffAnalysisRes)
%RADIALARRANGEMENT 'quantifies' the radiality of track arrangement and gives back the center
%
%SYNOPSIS [centerCoord,radialityParam] = radialArrangement(tracks,...
%    centerCoordInit,minTrackLen,diffAnalysisRes)
%
%INPUT  tracks         : Output of trackCloseGapsKalman.
%       centerCoordInit: Initial guess of center of radial track
%                        arrangement.
%                        Optional. Default: [0 0]
%       minTrackLen    : Minimum length of a track to be used in analysis.
%                        Optional. Default: 5.
%       diffAnalysisRes: Diffusion analysis results (output of
%                        trackDiffusionAnalysis1). Optional. If not input, will
%                        be calculated here.
%
%OUTPUT radialityParam: Parameter quantifying how radial the arrangment of
%                       tracks is. Calculated as the average angle (in
%                       radians) between the vector coming out of the "cell
%                       center" to each track center and the track direction.
%                       0 = perfect radiality. The less radial, the higher
%                       the value.
%       centerCoord   : Coordinates of center of radial track arrangement.
%       numTracksLin  : Number of linear tracks used in calculation.
%
%Khuloud Jaqaman, June 2009

%% input

if nargin < 2 || isempty(tracks) || isempty (centerCoordInit)
    centerCoordInit = [0 0];
end

if nargin < 3 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 4 || isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(tracks,1,2,1,[0.05 0.1]);
end

%% preamble

%ignore merges and splits and divide compound tracks back into the
%individual track segments
inputStruct = tracks;
clear tracksFinal
tracks = convStruct2MatIgnoreMS(inputStruct);

%extract track classification from diffAnalysisRes
trackType = vertcat(diffAnalysisRes.classification);

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
tracks = tracks(indx,:);
trackType = trackType(indx,:);

%find linear tracks
tracksLin = tracks(trackType(:,1)==1,:);
numTracksLin = size(tracksLin,1);

clear criteria tracks indx

%% calculation of center and radiality parameter

%find center coordinates such that the angles between rays coming out
%of the center and the track's direction of motion is minimized
options = optimset('MaxFunEvals',3000);
[centerCoord,radialityParam] = fminsearch(@angleRaysTracks,centerCoordInit',options,tracksLin);
centerCoord = centerCoord';


%% ~~~ the end ~~~


%% subfunction for radiality parameter calculation (and minimization)

function radialityParam = angleRaysTracks(centerCoord,tracks)

%get number of tracks
numTracks = size(tracks,1);

%reserve memory
dirAngle = NaN(numTracks,2);

%go over all tracks
for iTrack = 1 : numTracks

    %get the positions in this track and their standard deviations
    %keep NaNs to mark gaps
    trackCoordX = tracks(iTrack,1:8:end)';
    deltaCoordX = tracks(iTrack,5:8:end)';
    trackCoordY = tracks(iTrack,2:8:end)';
    deltaCoordY = tracks(iTrack,6:8:end)';
    trackCoord = [trackCoordX trackCoordY];
    deltaCoord = [deltaCoordX deltaCoordY];

    %estimate track's center
    trackCenter = nanmean(trackCoord);

    %estimate track's direction of motion
    [dummy,dummy,trackDir] = projectCoordOntoDir(trackCoord,deltaCoord);

    %calculate direction of vector from radial arrangement center to track
    %center
    rayDir = trackCenter - centerCoord';
    rayDir = rayDir / norm(rayDir);

    %calculate the angle between the two direction vectors
    tmp = min(abs(rayDir * trackDir),1);
    dirAngle(iTrack,1) = acos(tmp);
    
    %save also the number of points in this track
    dirAngle(iTrack,2) = length(find(~isnan(trackCoordX)));

end %(for iTrack = 1 : numTracks)

%calculate the weighted average of the angles (give longer tracks more
%weight)
radialityParam = sum(dirAngle(:,1).*dirAngle(:,2))/sum(dirAngle(:,2));

