function [data,label,legendText,selectedTags] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyz,legendText,correction)
%get tag number
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);

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
        
    case 1 %centroid correction
        %get centroids
        centroids = cat(1,anaDat.centroid);
        
        %correct coords
        tagCoord = tagCoord - centroids;
        
        correctionString = 'centroid corr.'
        
    case 2 %tag correction
        corrTag = get(handles.cenHandles(2),'Value')-1;
        
        corrTagCoord = squeeze(coordList(corrTag,:,:))';
        
        %correct doords
        tagCoord = tagCoord - corrTagCoord;
        
        correctionString = [anaDat(1).info.labelColor{corrTag},' centered'];
end

%calculate avg displacement for deltaT = 1...N/2: not fast to run,
%but fast to program
nDeltas = floor(timePoints(end)/2);
displacements = zeros(length(timePoints)-1,nDeltas);
deltas = displacements; %zeros, too

%calculate displacements
for deltaT = 1:nDeltas
    deltaCoord = tagCoord(1+deltaT:end,:) - tagCoord(1:end-deltaT,:);
    nDeltaCoord = length(deltaCoord);
    %store displacements.^2 and corresponding deltaT's
    displacements(1:nDeltaCoord,deltaT) = sum(deltaCoord.^2,2);
    deltas(1:nDeltaCoord,deltaT) = timePoints(1+deltaT:end) - timePoints(1:end-deltaT);
end

%calculate mean displacements corresponding to the timesteps stored
%in deltas
deltaVector = deltas(:);
deltaVector = unique(deltaVector);
deltaVector = deltaVector(2:end)';
%deltaVector now contains all timesteps that have been found
nDeltas = deltaVector(end);
meanDisplacement = zeros(nDeltas,1);

%read out mean
for deltaT = deltaVector
    meanDisplacement(deltaT) = mean(displacements(find(deltas == deltaT)));
end

data = meanDisplacement;

%write axis label
label = ['Displacement [\mum^2]'];

legendText = [legendText;{['Diffusion of ',anaDat(1).info.labelColor{tag1}, ' [microns^2] ',correctionString]}];