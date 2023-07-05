function [h] = overlayMaskBoundaries(mask,color)
%function [] = overlayMaskBoundaries(adhDilated) overlays the mask
%boundaries in the figure.
%   mask            mask
% Sangyoon Jul 4 2023

if nargin<2
    color = 'c';
end
curAdhBound = bwboundaries(mask,4);
for kk=1:numel(curAdhBound) 
    adhBoundary = curAdhBound{kk}; %(i).adhBoundary;
    h(kk) = line(adhBoundary(:,2), adhBoundary(:,1),...
                'Color',color,'LineStyle','-');
end

end

