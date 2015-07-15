function [] = UIPartitionAnalysis(track2mask)
% Determines the partitioing fraction of tracks in given MovieList
%
%First ML prompt needs ML with TrackingProcess. While, the second ML prompt
%needs that with DetectionProcess.
%
% SYNOPSIS function [] = UIPartitionAnalysis(track2mask)
%
% INPUT
%   track2mask      : indices of mask MovieList that corresponds to track
%                     MovieList. For example, MaskIndx = track2mask(TrackIndx),
%                     where MaskIndx corresponds to index of MD in mask ML
%                     and TrackIndx corresponds to index of MD in track ML.
%                     In other words, this is a map that translates track
%                     indices to mask indices.
%
%Tae H Kim, July 2015

%% User prompt
%ML for tracks
[fileNameTrack, filePathTrack] = uigetfile('*.mat', 'Select MovieList containing tracks for partition analysis');
%ML for masks
[fileNameMask, filePathMask] = uigetfile('*.mat', 'Select MovieList containing mask information for partition analysis');

%% Input Check
% the length of MLs and track2mask must be equal
%load the MLs
ML_Track = MovieList.load([filePathTrack fileNameTrack]);
ML_Mask = MovieList.load([filePathMask fileNameMask]);
%determines the length and compares them
nMD = numel(ML_Track.movies_);
if nMD ~= numel(ML_Mask.movies_)
    error('The size of ML containing track information does not match that of ML containing mask information.');
end
%assigns default value to track2mask
if isempty(track2mask)
    track2mask = 1 : nMD;
end
%compares the map length to ML
if nMD ~= numel(track2mask)
    error('The size of the map and the MovieLists do not match.')
end

%% Partition Analysis
for iMD = 1:nMD
    printLength = fprintf('MD %g/%g\n', iMD, nMD);
    maskDetectedStructure(ML_Mask.movies_{iMD});
    trackPartitioning(ML_Track.movies_{iMD}, ML_Mask.movies_{iMD});
    fprintf(repmat('\b', 1, printLength));
end

%% Save
ML_Mask.save();
ML_Track.save();

end

