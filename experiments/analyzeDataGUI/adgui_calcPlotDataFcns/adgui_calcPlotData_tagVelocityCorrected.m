function [data,label,legendText,selectedTags] = adgui_calcPlotData_tagVelocityCorrected(handles,anaDat,xyz,legendText,correction)
%get tag number
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);

selectedTags = [tag1];

%calc displacement vectors (see adgui_calcPlotData_tagDisplacementCorrected for comments)
displacementVectorsN = cat(3,anaDat.displacementVectorsN);
displacementVectorsN = squeeze(displacementVectorsN(tag1,:,:));
displacement = cat(2, anaDat.displacement);
displacement = ones(3,1) * displacement(tag1, :);
displacementVectors = displacementVectorsN .* displacement;

%switch on how to correct displacement
switch correction
    
    case 1 %centroid correction
        %get centroids
        centroids = cat(1,anaDat.centroid);
        
        %centroid displacement
        deltaCentroid = centroids(2:end,:)-centroids(1:end-1,:);
        
        %correct displacement
        correctedDisplacementVectors = displacementVectors' - deltaCentroid; 
        
        correctionString = '(centroid corr.)'
        
    case 2 %tag correction
        corrTag = get(handles.cenHandles(2),'Value')-1;
        
        coordList = cat(3,anaDat.coord);
        corrTagCoord = squeeze(coordList(corrTag,:,:))';
        
        deltaCorrTagCoord = corrTagCoord(2:end,:) - corrTagCoord(1:end-1,:);
        
        %correct displacement
        correctedDisplacementVectors = displacementVectors' - deltaCorrTagCoord; 
        
        correctionString = [anaDat(1).info.labelColor{corrTag},' centered'];
end

correctedDisplacement = normList(correctedDisplacementVectors);

time = cat(1,anaDat.time);
deltaT = (time(2:end) - time(1:end-1))/60;
velocityVector = correctedDisplacement./deltaT;
data = velocityVector;

%write axis label
label = ['Velocity [\mum/min]'];

legendText = [legendText;{['Velocity of ',anaDat(1).info.labelColor{tag1}, ' [microns/min] ',correctionString]}];
