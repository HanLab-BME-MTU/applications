function [data,label,legendText,selectedTags] = adgui_calcPlotData_tagDisplacementUncorrected(handles,anaDat,xyz,legendText)

%get tag number
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);

selectedTags = [tag1];

%find displacement
displacementMatrix = cat(2,anaDat.displacement);
data = displacementMatrix(tag1,:)';

%write axis label
label = ['Displacement [\mum]'];

legendText = [legendText;{['Displacement of ',anaDat(1).info.labelColor{tag1}, ' [microns] (uncorrected)']}];
