function [data,label,legendText,selectedTags] = adgui_calcPlotData_perpendicularProjection(handles,anaDat,xyz,legendText)
%v1 - v1//v2, where a//b is the projection of a onto b

%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);
eval(['tag3 = get(handles.',xyz,'TagHandles(3),''Value'')-1;']);
eval(['tag4 = get(handles.',xyz,'TagHandles(4),''Value'')-1;']);

selectedTags = [tag1,tag2,tag3,tag4];

%find vectors
distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
vector1 = squeeze(distanceVectorMatrixN(tag1,tag2,:,:))';
vector2 = squeeze(distanceVectorMatrixN(tag3,tag4,:,:))';

%find vector length
distanceMatrix = cat(3,anaDat.distanceMatrix);
vector1Length = squeeze(distanceMatrix(tag1,tag2,:));

%calculate projection
projVector = ((vector1Length .* dot(vector1,vector2,2)) * ones(1,3)) .* vector2;

%calculate perpendicular vector
perpVector = vector1.*(vector1Length*ones(1,3)) - projVector;

%calculate length
perpDistance = normList(perpVector);


data = perpDistance;

%write axis label
label = ['Distance [\mum]'];

legendText = [legendText;{['perp. dist. of proj. of',anaDat(1).info.labelColor{tag1}, '-', anaDat(1).info.labelColor{tag2},...
                ' onto ', anaDat(1).info.labelColor{tag3}, '-', anaDat(1).info.labelColor{tag4}, '[microns]']}];
