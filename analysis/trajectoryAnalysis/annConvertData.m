function taInput = annConvertData(annData,micronsPerPixel, uncertainty)
%ANNCONVERTDATA converts Ann's data so that trajectoryAnalysis can use it as first input argument
%
% SYNOPSIS: taInput = annConvertData(annData, micronsPerPixel)
%
% INPUT annData: n-by-4 array (further columns optional, but not considered)
%           C1: MT number
%           C2: x (pix)
%           C3: y (pix)
%           C4: sampling intervall (s) 
%       micronsPerPixel : pixelsize in microns
%       uncertainty (opt): uncertainty of the position measurements in
%                          pixels. Default: 1
%
% OUTPUT taInput: data structure to be used as first input argument in
%        trajectoryAnalysis 
%        taInput(1:n) : structure containing n different trajectories with fields
%           - distance   tx2 array [distance, sigmaDistance] in microns
%           - time       tx2 array [time, sigmaTime] in seconds
%           - timePoints tx1 array [timePoint#] 
%           - info       struct containing additional information about the inputData
%                  -tags    : cell containing the two strings designating
%                             the tags between which the distance is
%                             measured, e.g. {'spb1','cen1'}
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 14-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===========================
%% TEST INPUT
%===========================

% default: uncertainty of 1 pixel
def_uncertainty = 1;

if nargin < 2 || isempty(annData) || isempty(micronsPerPixel)
    error('annConvertData requires two non-empty input arguments')
end

if nargin < 3 || isempty(uncertainty);
    uncertainty = def_uncertainty;
end

%===========================


%===========================
%% CONVERT
%===========================

% find number of MTs 
nMTs = max(annData(:,1));

% preassign output

taInput(1:nMTs) = struct('distance',[],'time',[],...
    'timePoints',[],'info',struct('tags','singleMT'));

% loop through MTs and calculate distance and uncertainty. Since we have
% 2-d data, we use the dot product of consecutive displacement vectors to
% test whether we're continuing growth in the same direction, or whether
% there is a switch.

for iMT = 1:nMTs
    % find rows describing current MT
    currentMT = annData(:,1)==iMT;
    nTimepoints = sum(currentMT);
    
      
    % read xy-positions. Convert to microns. Assume that pixels are square.
    points.coordinates = annData(currentMT,2:3) .* micronsPerPixel;
    % covariances. Meaningless at this point. We will add uncertainty to
    % the eventual MT length.
    points.covariances = ones(size(points.coordinates));
    
    % calculate displacement
    [displacements, dummy, unitVector] = deltaCoordinates(points);
    
    % assign direction, using scalar product of unitVectors
    direction = sum(unitVector(2:end,:).*unitVector(1:end-1,:),2);
    direction = sign(direction);
    % !!!! assign growth to sign 0 
    direction(direction == 0) = 1;
    displacements = displacements .* [1;direction];
    
    % turn displacements into distance measurements
    distance = cumsum([0;displacements]);
    
    % shift the distance so that we don't go below 0
    % distance = distance - min(distance) +
    % 0.1*(max(distance)-min(distance))
    distance = distance - 1.1*min(distance) + 0.1*max(distance);
    
    % save distance, add uncertainty
    taInput(iMT).distance = ...
        [distance,uncertainty*micronsPerPixel*ones(size(distance))];
    
    % add timepoints. Assume for the moment that there are no lost frames
    taInput(iMT).timePoints = (1:nTimepoints)';
    
    % add time. Ann puts 0 at the end of each MT-set
    time = annData(currentMT,4);
    taInput(iMT).time = [0;cumsum(time(1:end-1))];
end