function cleanFocalAdhesionMaskSample = cleanUpFocalAdhesionMasks( focalAdhesionMaskSample )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

focalAdhesionMaskSample=(focalAdhesionMaskSample > 0);
labels=bwlabel(focalAdhesionMaskSample);
stats = regionprops(labels, 'Orientation','Area');
idx = find(([stats.Area] > 10)&(abs([stats.Orientation])<20));
cleanFocalAdhesionMaskSample=ismember(labels, idx);

end

