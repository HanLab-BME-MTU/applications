function [data,label,legendText,selectedTags,stats] = adgui_calcPlotData_distance(handles,anaDat,xyz,legendText,dataProperties)

%distance between tags
%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);

selectedTags = [tag1,tag2];

%find distances
distanceMatrix=cat(3,anaDat.distanceMatrix);
data = squeeze(distanceMatrix(tag1,tag2,:));

%find distance vectors (to calculate sigma)
distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
distanceVectors = squeeze(distanceVectorMatrixN(tag1,tag2,:,:))' .* (data*ones(1,3));

%write axis label
labelString = 'Distance [\mum]';
label = labelString;

legendText = [legendText;{['distance ', anaDat(1).info.labelColor{tag1},...
                '-', anaDat(1).info.labelColor{tag2}, ' [microns]']}];

[stats.sigma] = adgui_calcPlotData_distanceSigma(anaDat,selectedTags,distanceVectors,data,0);