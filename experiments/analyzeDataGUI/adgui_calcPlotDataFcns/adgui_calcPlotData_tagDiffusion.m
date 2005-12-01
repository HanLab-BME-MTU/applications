function [data,label,legendText,selectedTags,stats,time] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyz,legendText,correction,dataProperties,verbose)
%calculate diffusion data in the analyzeData framework. 
%
% SYNOPSIS [data,label,legendText,selectedTags,stats] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyz,legendText,correction)
%
% INPUT    handles:     either handles to analyzeDataGUI, or vector, where
%                       handles (1) is the number of the first tag and handles(2) is the
%                       correction tag (if there's one)
%	       anaDat:      basic data structure, created with calcAnaDat
%	       xyz:         string ('x','y','z'), depending on which axis to plot in (can be empty)
%	       lengendText: cell text array, containing legend text (can be empty)
%          correction:  how should tag positions be corrected? 
%		                0: no correction
%		                1: centroid correction
%		                2: tag correction
%                       3: 1D-diffusion
%          dataProperties: needed if you want to display time, else
%                           optional
%          verbose      whether or not to show popups (default 1, for analyzeDataGUI)
%
% OUTPUT    data:           mean squared diffused distance (during 1,2,3 etc. timesteps)
%	        label:          data label
%	        legendText:     legend text
%	        selectedTags:   number of diffusing tag
%	        stats:          structure, stats.sigma gives precision for each measurement
%           time:           time associated with the respective deltas. if
%                           this is returned, the code only gives
%                           data, sigma, where there are no NaNs!
%
%c: 6/03, jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargout > 5 & nargin < 6 | isempty(dataProperties)
    error('Please specify dataProperties')
end

if nargin < 7 | isempty(verbose)
    verbose = 1;
end

%get tag number
if isstruct(handles)
    eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
else
    tag1 = handles(1);
end

selectedTags = [tag1];

%get timepoints
timePoints = cat(1,anaDat.timePoint);

%get coords
coordList = cat(3,anaDat.coord);
tagCoord = squeeze(coordList(tag1,:,:))'; %tagCoord is a Nx3 matrix

%switch on how to correct diffusion
switch correction
    
    case 0 %no correction
        correctionString = 'uncorrected';
        
        %set indicate that there is no correction
        corrTag = 0;
        
    case 1 %centroid correction
        %get centroids
        centroids = cat(1,anaDat.centroid);
        
        %correct coords
        tagCoord = tagCoord - centroids;
        
        correctionString = 'centroid corr.'
        
        %set the correcting number to nTags+1 (which is where the q-matrix is stored)
        corrTag = anaDat(1).info.nTags + 1;
        
    case 2 %tag correction
        
        if isstruct(handles)
            corrTag = get(handles.cenHandles(2),'Value')-1;
        else
            corrTag = handles(2);
        end
        
        corrTagCoord = squeeze(coordList(corrTag,:,:))';
        
        %correct coords
        tagCoord = tagCoord - corrTagCoord;
        
        correctionString = [anaDat(1).info.labelColor{corrTag},' centered'];
        
    case 3 %1-D diffusion, tag corrected
        
        if isstruct(handles)
            corrTag = get(handles.cenHandles(2),'Value')-1;
        else
            corrTag = handles(2);
        end
        
        % new tagCoord will be: [distance, 0, 0], which will effectively
        % turn this into a 1D-diffusion
        
        %get distance
        distanceMatrix=cat(3,anaDat.distanceMatrix);
        distance = squeeze(distanceMatrix(tag1,corrTag,:));
        
        %correct coords
        tagCoord = zeros(size(distance,1),3);
        tagCoord(:,1) = distance;
        
        correctionString = [anaDat(1).info.labelColor{corrTag},' centered (1-D)'];
end

%calculate avg displacement for deltaT = 1...N/4: not so fast to run,
%but fast to program
%do not use more than 1/4 of the data. See Saxton 1997, Biophys. J 
%use length(timepoints), because there could be several tps missing!
nDeltas = floor(length(timePoints)/4);
%check that we have enough data
if nDeltas<1
    if verbose
        h = errordlg('not enough data points for diffusion measurement');
        uiwait(h)
        return
    else
        error('not enough data points for diffusion measurement');
    end
end
displacements = zeros(length(timePoints)-1,nDeltas);
displacementSigma = displacements; %init sigmas the same way, so that both matrices have the same size!
deltas = displacements; %zeros, too

%calculate displacements. We have to calculate sigmas here (not sorted by deltas yet)
%because deltaT is measured in entries in anaDat and not in real time 

for deltaT = 1:nDeltas
    %calculate displacements - then check for 0-displacements (can occur for fusions!)
    displacementVectors = tagCoord(1+deltaT:end,:) - tagCoord(1:end-deltaT,:);
    displacementNorms2 = sum(displacementVectors.^2,2);
    badVectors = find(displacementNorms2==0);
    
    %calculate sigmas first with all vectors, then discard bad data
    nDisplacementVectorsAll = length(displacementNorms2);
    
    %calculate sigmas
    displacementSigmaTmp =...
        adgui_calcPlotData_distanceSigma(anaDat,selectedTags,displacementVectors,sqrt(displacementNorms2),corrTag,deltaT);
    
    %throw away bad data
    displacementSigmaTmp(badVectors) = [];
    displacementNorms(badVectors) = [];
    displacementVectors(badVectors,:) = [];
    
    %get number of good vectors
    nDisplacementVectors = length(displacementNorms2);
    
    %store displacements.^2 
    displacements(1:nDisplacementVectors,deltaT) = displacementNorms2;
    
    %store sigmas
    displacementSigma(1:nDisplacementVectors,deltaT) = displacementSigmaTmp;
    
    %store deltas corresponding to the displacements. Account for discarded
    %displacements
    deltaTmp = timePoints(1+deltaT:end) - timePoints(1:end-deltaT);
    deltaTmp(badVectors) = [];
    deltas(1:nDisplacementVectors,deltaT) = deltaTmp;
end

%calculate mean displacements corresponding to the timesteps stored
%in deltas
deltaVector = deltas(:);
deltaVector = unique(deltaVector);
deltaVector = deltaVector(2:end)'; %use 2:end because deltaVector contains zeros
%deltaVector now contains all timesteps that have been found
nDeltas = deltaVector(end);

%----read out statistics
%init variables: fill with NaNs
[meanDisplacement,sigmaDisplacement] = deal(repmat(NaN,nDeltas,1));

for deltaT = deltaVector
    deltaIdx = find(deltas == deltaT);
    %make sure we use enough data (since there might have been deleted
    %frames, the last entries might not bee good enough)
    if length(deltaIdx)>20
        [meanDisplacement(deltaT),sigmaDisplacement(deltaT)] = weightedStats(displacements(deltaIdx),displacementSigma(deltaIdx),'s');
    else
        [meanDisplacement(deltaT),sigmaDisplacement(deltaT)] = deal(NaN);
    end
end

%test whether there is at least 1 good timepoint
if all(isnan(meanDisplacement+sigmaDisplacement))
    if verbose
        h = errordlg('there is not one valid data point (we need at least 20 measurements for one data point)');
        uiwait(h)
    end        
    %return does not work here, so we rethrow the error (why??)
    error('no data to plot');
end

data = meanDisplacement;
stats.sigma = sigmaDisplacement;

%write axis label
label = ['Displacement [\mum^2]'];

legendText = [legendText;{['Diffusion of ',anaDat(1).info.labelColor{tag1}, ' [microns^2] ',correctionString]}];

%---return timesteps if asked for
if nargout > 5
    % timepoints2secondsFactor
    timeFactor = dataProperties.timeLapse;
    
    % only return time, data, where there are actually values
    timePoints = [1:length(data)];
    nanIdx = find(isnan(data));
    data(nanIdx) = []; %kill nans
    stats.sigma(nanIdx) = []; %kill nans
    timePoints(nanIdx) = [];
    
    time = timeFactor*timePoints;
end