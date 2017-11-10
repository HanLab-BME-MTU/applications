function bw_out = imdilateWithScale(bw_in,scaleMap,levels)
% dilate a bw image with the corresponding scale

bw_out = bw_in;

for i = 1 : length(levels)
    if levels(i)>0
        H = fspecial('disk', ceil(levels(i)))>0;
        this_scale_bw = and(scaleMap==i, bw_in>0);
        this_scale_bw = imdilate(this_scale_bw, H,'same');
        bw_out = or(bw_out,this_scale_bw);
    end
end
