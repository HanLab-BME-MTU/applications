function [data,label,legendText,selectedTags] = adgui_calcPlotData_angleVectors(handles,anaDat,xyz,legendText)
%angle between vectors

%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);
eval(['tag3 = get(handles.',xyz,'TagHandles(3),''Value'')-1;']);
eval(['tag4 = get(handles.',xyz,'TagHandles(4),''Value'')-1;']);

selectedTags = [tag1,tag2,tag3,tag4];

%find vectors
distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
vector1 = squeeze(distanceVectorMatrixN(tag1,tag2,:,:));
vector2 = squeeze(distanceVectorMatrixN(tag3,tag4,:,:));

angles = acos(sum(vector1.*vector2,1)')*180/pi;
data = angles;

%write axis label
label = ['angle [deg]'];

legendText = [legendText;{['angle between ',anaDat(1).info.labelColor{tag1}, '-', anaDat(1).info.labelColor{tag2},...
                ' and ', anaDat(1).info.labelColor{tag3}, '-', anaDat(1).info.labelColor{tag4}, '[rad]']}];