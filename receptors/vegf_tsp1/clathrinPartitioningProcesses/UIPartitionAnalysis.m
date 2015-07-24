function [] = UIPartitionAnalysis(varargin)
% Determines the partitioing fraction of tracks in given MovieList
%
%First ML prompt needs ML with TrackingProcess. While, the second ML prompt
%needs that with DetectionProcess.
%
% SYNOPSIS function [] = UIPartitionAnalysis(track2mask)
%
% INPUT
%   varargin        : name value pair
%       'psfSigmaMult'          : ultiplicative factor used to determine the size
%                                 of mask. radius = psfSigmaMult *
%                                 pasfSigma. Part of maskDetectedStructure.
%       'scrambleTracks'        : scrambles mean track position for control
%                                 data set generation
%
%Tae H Kim, July 2015

%% Input
%input parser
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('psfSigmaMult', 3, @isnumeric);
ip.addParameter('scrambleTracks', false, @(x) islogical(x) || isnumeric(x));
ip.parse(varargin{:});
psfSigmaMult = ip.Results.psfSigmaMult;
scrambleTracks = ip.Results.scrambleTracks;
fprintf('Tracks will be scrambled');

%% User prompt
%ML for tracks
[fileNameTrack, filePathTrack] = uigetfile('*.mat', 'Select MovieList containing tracks for partition analysis'); %#ok<ASGLU>
%ML for masks
[fileNameMask, filePathMask] = uigetfile('*.mat', 'Select MovieList containing mask information for partition analysis'); %#ok<ASGLU>

%% Input Check
% the length of MLs and track2mask must be equal
%load the MLs
pLength = fprintf('Loading ML_track\n');
evalc('ML_Track = MovieList.load([filePathTrack fileNameTrack])');
fprintf(repmat('\b', 1, pLength));
fprintf('Loading ML_track complete\n');
pLength = fprintf('Loading ML_mask\n');
evalc('ML_Mask = MovieList.load([filePathMask fileNameMask])');
fprintf(repmat('\b', 1, pLength));
fprintf('Loading ML_mask complete\n');
%determines the length and compares them
nMD = numel(ML_Track.movies_);
if nMD ~= numel(ML_Mask.movies_)
    error('The size of ML containing track information does not match that of ML containing mask information.');
end
%{
%assigns default value to track2mask
if nargin<1
    track2mask = 1 : nMD;
end
%compares the map length to ML
if nMD ~= numel(track2mask)
    error('The size of the map and the MovieLists do not match.')
end
%}

%% Partition Analysis
multipleProgressText('Analyzing MD', nMD);
for iMD = 1:nMD
    %SubcellMaskProcess-------------------------------------------------------
    %get default para
    maskPara = SubcellMaskProcess.getDefaultParams(ML_Mask.movies_{iMD}.outputDirectory_);
    maskPara.psfSigmaMult = psfSigmaMult;
    %call analysis function
    maskDetectedStructure(ML_Mask.movies_{iMD}, maskPara);
    %track partitioning process--------------------------------------------
    %get default para
    trackPara = PartitionAnalysisProcess.getDefaultParams(ML_Mask.movies_{iMD}.outputDirectory_);
    trackPara.scrambleTracks = scrambleTracks;
    %call analysis function
    trackPartitioning(ML_Track.movies_{iMD}, ML_Mask.movies_{iMD}, trackPara);
    multipleProgressText();
end

%% Save
ML_Mask.save();
ML_Track.save();

end

