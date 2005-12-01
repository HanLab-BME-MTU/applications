function [data,label,legendText,selectedTags] = adgui_calcPlotData_deltaDistanceVelocity(handles,anaDat,xyz,legendText)

%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);

selectedTags = [tag1,tag2];

%find distances
distanceMatrix=cat(3,anaDat.distanceMatrix);
distance = squeeze(distanceMatrix(tag1,tag2,:));
deltaDistance = distance(2:end) - distance(1:end-1);

time = cat(1,anaDat.time);
deltaT = (time(2:end) - time(1:end-1))/60;
velocityVector = abs(deltaDistance./deltaT);
data = velocityVector;

%write axis label
label = 'velocity [\mum/min]';

legendText = [legendText;{['delta distance(', anaDat(1).info.labelColor{tag1},...
                '&', anaDat(1).info.labelColor{tag2}, ')/time  [microns/min]']}];