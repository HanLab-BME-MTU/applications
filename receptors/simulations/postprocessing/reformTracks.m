function [ tracksNew] = reformTracks( compTracks,reformParam )
%REFORMTRACKS this function reform the compTracks for the simulated data,
%to have the same size and time step as the tracks from U-track analysis.
%SYNOPSIS tracksNew = reformTracks( compTracks,reformParam )
% Inputs
%       compTracks      : Tracks from the simulations.
%      
%       reformParam     : Structure with fields:
%                       
%          .simTimeStep     : Time step in simulated tracks.
%          .sampleTimeStep  : Sampling time step, in same units as timeStep.
%          .frameRange      : Row vector with two entries indicating time range
%                             to retain.
%           .pixelSize      : pixel size in microns
%           .psfSigma       : width of point-spread function (in pixels)
%
%
% 
%OUTPUT
%      tracksNew: reformed compTracks
%       
%
% Luciana de Oliveira, October 2016.
%% Input

%call imput variables from structure

% time info
simTimeStep = reformParam.simTimeStep;
sampleTimeStep = reformParam.sampleTimeStep;
frameRange = reformParam.frameRange;

%space info
pixelSize =reformParam.pixelSize;
psfSigma =reformParam.psfSigma;

%% Reformat ground truth tracks for comparison
tracks0 = sparseCompTracks(compTracks);
tracksSub = subSampleCompTracks(tracks0,simTimeStep,sampleTimeStep);
tracksSubCrop = cropCompTracksFR(tracksSub,frameRange);
tracksSim = tracksSubCrop;

%rescaling in function of pixel size
for iTrack = 1 : length(tracksSim)
    tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end) = tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end)/pixelSize;
    tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end) = tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end)/pixelSize; 
end
%changing the coordinates of tracksSim
for iTrack = 1 : length(tracksSim)
    tmp = tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end);
    tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end) = tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end);
    tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end) = tmp;
end

% Alleviate boundary effects on detection
for iTrack = 1 : length(tracksSim)
    tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end) = tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end)+10*psfSigma;
    tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end) = tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end)+10*psfSigma;
end
tracksNew = sepCompTracks(tracksSim);


end

