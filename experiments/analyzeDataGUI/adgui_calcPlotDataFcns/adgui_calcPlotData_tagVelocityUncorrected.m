function [data,label,legendText,selectedTags] = adgui_calcPlotData_tagVelocityUncorrected(handles,anaDat,xyz,legendText)
%calculate tag velocity

%get tag number
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);

selectedTags = [tag1];

%calc velocity
displacementMatrix = cat(2,anaDat.displacement);
tagDisplacement = displacementMatrix(tag1,:)';
time = cat(1,anaDat.time);
deltaT = (time(2:end) - time(1:end-1))/60;
velocityVector = tagDisplacement./deltaT;
data = velocityVector;

%write axis label
label = ['Velocity [\mum/min]'];

legendText = [legendText;{['Velocity of ',anaDat(1).info.labelColor{tag1}, ' [microns/min] (uncorrected)']}];