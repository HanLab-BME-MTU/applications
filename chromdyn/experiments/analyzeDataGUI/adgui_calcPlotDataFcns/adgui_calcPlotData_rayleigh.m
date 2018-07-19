function [data,label,legendText,selectedTags] = adgui_calcPlotData_rayleigh(handles,anaDat,xyz,legendText,dataProperties)

%rayleigh (distance)

%get tag numbers
eval(['tag1 = get(handles.',xyz,'TagHandles(1),''Value'')-1;']);
eval(['tag2 = get(handles.',xyz,'TagHandles(2),''Value'')-1;']);

selectedTags = [tag1,tag2];

%find vector
distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
vector1 = squeeze(distanceVectorMatrixN(tag1,tag2,:,:));

sc = dataProperties.sigmaCorrection;

%find rayleigh-distances for e1,e2,e3
rayleighVector = dataProperties.FT_SIGMA.*[0.61/0.21/sc(1), 0.61/0.21/sc(1), 2/0.66/sc(2)];
%calc in microns
rayleighVector = rayleighVector .* [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];

%calculate rayleigh-distance: radius of the ellipsoid defined by the
%rayleigh-distances in the main axis in the direction of vector1
rayleighDistances = ellipsoidRadius(vector1',rayleighVector);

data = rayleighDistances'';

%write axis label
labelString = 'Distance [\mum]';
label = labelString;

legendText = [legendText;{['rayleigh-limit for vector ', anaDat(1).info.labelColor{tag1},...
                '-', anaDat(1).info.labelColor{tag2}, ' [microns]']}];