function out_img = keep_area_with_threshold(in_img,T_low, T_high)

labelMask = bwlabeln(in_img);

if ~exist('T_high','var')
    T_high = numel(in_img(:));
end

%Get their area
obProp = regionprops(labelMask,'Area','BoundingBox');
out_img=zeros(size(in_img))>0;

%First, check that there are objects to remove
if length(obProp) > 1
    obAreas = [obProp.Area];
    obBoundingBox = reshape([obProp.BoundingBox],[6 numel(obAreas)]);
    zWidth = obBoundingBox(6,:);
    ind = find(obAreas>T_low & obAreas<=T_high & zWidth>2 );
    
    for i = ind
        out_img = or(out_img,labelMask == i);
    end
else
    out_img = in_img;
end
        