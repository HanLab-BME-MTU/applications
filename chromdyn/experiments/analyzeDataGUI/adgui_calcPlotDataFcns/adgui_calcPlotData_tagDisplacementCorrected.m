function [data,label,legendText,selectedTags,stats] = adgui_calcPlotData_tagDisplacementCorrected(handles,anaDat,xyz,legendText,dataProperties,correction,tagNumbers)

%check input: if no handles passed down, we use the tagNumbers
tagsFromInput = 0;
if isempty(handles)
    if nargin > 6 & ~isempty(tagNumbers)
        tag1 = tagNumbers(1);
        selectedTags = tag1;
        tagsFromInput = 1;
    else
        error('not enough input arguments: either handles or tagNumber are missing!');
    end
else
    %get tag number
    eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
    
    selectedTags = [tag1];
end

%find displacementVectors
displacementVectorsN = cat(3,anaDat.displacementVectorsN);

%choose the right one
displacementVectorsN = squeeze(displacementVectorsN(tag1,:,:));

%find vector lengths
displacement = cat(2, anaDat.displacement);

%prepare matrix for multiplication
displacement = ones(3,1) * displacement(tag1, :);

%reconstruct original displacementVectors
displacementVectors = displacementVectorsN .* displacement;

%switch on how to correct displacement
switch correction
    
    case 0 %no correction
        
        %no calcs, just write out results
        correctedDisplacementVectors = displacementVectors';
        correctionString = '(uncorrected)';
        
        correctionTag = 0;
        
    case 1 %centroid correction
        %get centroids
        centroids = cat(1,anaDat.centroid);
        
        %centroid displacement
        deltaCentroid = centroids(2:end,:)-centroids(1:end-1,:);
        
        %correct displacement
        correctedDisplacementVectors = displacementVectors' - deltaCentroid; 
        
        correctionString = '(centroid corr.)';
        
        correctionTag = anaDat(1).info.nTags+1;
        
    case 2 %tag correction
        if tagsFromInput
            correctionTag = tagNumbers(2);
        else
            correctionTag = get(handles.cenHandles(2),'Value')-1;
        end
        
        coordList = cat(3,anaDat.coord);
        corrTagCoord = squeeze(coordList(correctionTag,:,:))';
        
        deltaCorrTagCoord = corrTagCoord(2:end,:) - corrTagCoord(1:end-1,:);
        
        %correct displacement
        correctedDisplacementVectors = displacementVectors' - deltaCorrTagCoord; 
        
        correctionString = [anaDat(1).info.labelColor{correctionTag},' centered'];
        
end


%calc norms
correctedDisplacement = normList(correctedDisplacementVectors);

data = correctedDisplacement;

%write axis label
label = ['Displacement [\mum]'];

legendText = [legendText;{['Displacement of ',anaDat(1).info.labelColor{tag1}, ' [microns] ',correctionString]}];

[stats.sigma] = adgui_calcPlotData_distanceSigma(anaDat,selectedTags,correctedDisplacementVectors,correctedDisplacement,correctionTag);