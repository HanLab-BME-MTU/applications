function [data,label,legendText,selectedTags,sigma] = adgui_calcPlotData_deltaDistance(handles,anaDat,xyz,legendText)
%delta distance

%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);

selectedTags = [tag1,tag2];

%find distances
distanceMatrix=cat(3,anaDat.distanceMatrix);
distance = squeeze(distanceMatrix(tag1,tag2,:));

deltaDistance = distance(2:end) - distance(1:end-1);

data = deltaDistance;

%write axis label
label = 'Distance [\mum]';

legendText = [legendText;{['delta distance ', anaDat(1).info.labelColor{tag1},...
                '-', anaDat(1).info.labelColor{tag2}, ' [microns]']}];

%extract data which is to be passed down for the sigma calculation
distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
distanceVector = squeeze(distanceVectorMatrixN(tag1,tag2,:,:))' .* (distance*ones(1,3));
%calculate sigma
[sigma] = adgui_calcPlotData_deltaDistanceSigma(anaDat,tags,distanceVectors,distances);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

