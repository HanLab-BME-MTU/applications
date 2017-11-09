function  new_OtsuRosin_Segment = ...
    nucleus_segmentation_refinemask(OtsuRosin_Segment,min_size, max_num, close_radius)

% close first
seClose = strel('disk',close_radius,0);
OtsuRosin_Segment = imclose(OtsuRosin_Segment, seClose);

% fill hole without the edge effect
OtsuRosin_Segment_pad = zeros(size(OtsuRosin_Segment)+2);
OtsuRosin_Segment_pad(2:end-1,2:end-1)=OtsuRosin_Segment;

OtsuRosin_Segment_pad = imfill(OtsuRosin_Segment_pad,'hole');
OtsuRosin_Segment = OtsuRosin_Segment_pad(2:end-1,2:end-1);


% for keeping which ones
labelMask= bwlabel(OtsuRosin_Segment);

obAreas = regionprops(labelMask,'Area');       %#ok<MRPBW>
obAreas = [obAreas.Area];

%First, check that there are objects to remove
if length(obAreas) > max_num || min(obAreas)<min_size
    [dummy,iSort] = sort(obAreas,'descend'); %#ok<ASGLU>
    %Keep only the largest requested number
    OtsuRosin_Segment = false(size(OtsuRosin_Segment));
    for i = 1:min(max_num,length(obAreas) )
        if(obAreas(iSort(i))>min_size)
            OtsuRosin_Segment = OtsuRosin_Segment | labelMask == iSort(i);
        end
    end
end

% separate touching ones

invert_seg = 1-OtsuRosin_Segment;
D = bwdist(invert_seg);
D = imfilter(D, fspecial('gaussian',100,15),'same','replicate');
invert_D = -D;
L = watershed(invert_D);
new_OtsuRosin_Segment = (double(L).*double(OtsuRosin_Segment))>0;



      