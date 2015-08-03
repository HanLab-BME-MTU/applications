function out_img = keep_largest_area(in_img)

labelMask = bwlabel(in_img);

%Get their area
obAreas = regionprops(labelMask,'Area');
out_img=in_img;

%First, check that there are objects to remove
if length(obAreas) > 1
    obAreas = [obAreas.Area];
    %Sort by area
    [dummy,iSort] = sort(obAreas,'descend');
    %Keep only the largest requested number
    out_img = labelMask == iSort(1);
else
    out_img = in_img;
end
        