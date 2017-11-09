function bw_out = gcaImdilateWithScale(bw_in,scaleMap,levels)
% dilate a bw image with the corresponding scale

bw_out = bw_in;
bw_in = logical(bw_in); 

%% SMALL quick and dirty to not trust the largest scale size
% get the highest level nearest neighbor and set it to that value  
% pix_high = find((scaleMap == levels(end) | scaleMap == levels(end-1))& bw_in>0) ; 
% bw_inNoHi  = bw_in; 
% bw_inNoHi(pix_high) = 0; 
% scaleMed = median(scaleMap(bw_inNoHi)); 
% if ~isinteger(scaleMed/0.5); 
%     scaleMed = round(scaleMed); 
% end 
% if isnan(scaleMed); 
%     % guess size 
%     scaleMed = 3; 
% end 
% scaleMap(pix_high) = scaleMed; 

% ideally want the pixels right before the pixel turns to noise 

%%

for i = 1 : length(levels)-2 % don't trust max level
    if levels(i)>0
     % H = fspecial('gaussian',[9,9], levels(i))>0;
       H = fspecial('disk',levels(i))>0; 
        this_scale_bw = and(scaleMap==levels(i), bw_in>0);
        this_scale_bw = imdilate(this_scale_bw, H,'same');
        bw_out = or(bw_out,this_scale_bw);
    end
end
