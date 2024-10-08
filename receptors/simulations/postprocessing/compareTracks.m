% function [ meanLinkStatsPerc,meanGapStatsPerc,meanMergeStats,meanSplitStats,mergeSplitStats] = compareTracks(compTracks,tracksFinal,reformParam)
function [linkStats,gapStats,meanLinkStatsPerc,meanGapStatsPerc] = compareTracks(compTracks,tracksFinal,reformParam)

%COMPARETRACKS this function compare the tracks from simulated data and
%uTrackPackageGUI.
%
%SYNOPSIS tracksOut = sparseCompTracks(tracksIn)
% Inputs
%       compTracks      : Tracks from the simulations.
%       tracksFinal     : Tracks from uTrackPackageGUI, in form of output of
%                         trackCloseGapsKalman.
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
%OUTPUT (as function SCORELINKSGAPSMS)
%       linkStats   : (Number of frames - 1) - by - 5 array. Row i
%                     corresponds to the links from frame i to frame i+1.
%                     The columns show the number of:
%                     (1) ground truth links;
%                     (2) links resulting from tracking;
%                     (3) tracking links that are correct;
%                     (4) tracking links that are between real features but
%                         are wrong;
%                     (5) tracking links that involve a detection artifact
%                         (i.e. a detetion false positive).
%                     The sum of columns 3 through 5 = column 2.
%       gapStats    : Two-column array with number of rows = number of
%                     gaps. First column shows gap length. Second column is
%                     0 - if gap is correctly closed;
%                     1 - if features just before and just after the gap
%                     are real but the connection is wrong;
%                     2 - if one or both features just before and just
%                     after the gap are detection artifacts.
%       mergeSplitStats: 2-by-3 matrix where first row is for merges, 2nd
%                        row is for splits.  1st column shows number of
%                        merges/splits in ground truth, 2nd column shows
%                        number of merges/splits in tracking results, and
%                        3rd column shows number of correct merges/splits
%                        in tracking results. Merge/split statistics
%                        don't explicitly account for detection false
%                        positives, i.e. they simply count toward wrong
%                        merges/splits, just like real features if wrongly
%                        assigned a merge/split.
%
% Luciana de Oliveira, September 2016.
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
%% Reformat U-tracks for comparison

%% Comparison 
[linkStats,gapStats] = scoreLinksGapsMS(tracksFinal,tracksNew);
meanLinkStats=mean(linkStats(:,3));
meanLinkStatsPerc=(100*meanLinkStats)/mean(linkStats(:,2));
totalNumRigthGaps=gapStats(gapStats(:,2)==0);
meanGapStatsPerc=(100*length(totalNumRigthGaps))/length(gapStats(:,2));
% meanMergeStats=mergeSplitStats(1,3);
% meanSplitStats=mergeSplitStats(2,3);
end %function

