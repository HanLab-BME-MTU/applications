function [mask, isClosed, outline, isClockWise] = getOutlineInfo(mask, iFrame)
%Load and combine masks from all channels
isClosed = ~any([mask(1,:), mask(end,:) mask(:,1)' mask(:,end)']);

%Separate mask objects - prSamProtrusion only handles one object
CC = bwconncomp(mask);

%Make sure there isn't more than one object
if CC.NumObjects > 1
%     if nChan == 1
        %Because small objects can be created by intersecting masks, we
        %only warn the user if only one mask was used.
        warning('blackwindow:protrusion:TooManyObjects',...
            ['The mask for frame ' num2str(iFrame) ...
            ' contains more than 1 object - protrusion vectors are only calculated for the largest object.'])
%     end
    %We just take the largest of these objects for prot vec calc
    [~,iBiggest] = max(cellfun(@(x)(numel(x)),CC.PixelIdxList));
    mask = false(size(mask));
    mask(CC.PixelIdxList{iBiggest}) = true;
    
end

%Get the outline of the object in this mask. We use contourc instead of
%bwboundaries for 2 reasons: It returns its results in matrix
%coordinates, and the resulting outline encloses the border pixels
%instead of running through their centers. This better agrees with the
%windows, as the windows are designed to enclose the entire mask.
mask(1,:)=0;
outline = contourc(double(mask),[0 0]);
outline = separateContours(outline);%Post-processing of contourc output
outline = cleanUpContours(outline);
outline = outline{1}';%We know we only have one object...

%Make sure the outline is correctly oriented
if ~isClosed
    %Close the curve before checking handedness
    closedOutline = closeContours({outline'},bwdist(~mask));
    isClockWise = isCurveClockwise(closedOutline{1});
else
    isClockWise = isCurveClockwise(outline);
end

if ~isClockWise
    %Sam requires the curves run in the same direction
    outline = outline(end:-1:1,:);
end
end