function [ imMultichannelOverlay ] = genMultiChannelOverlay( im, displaycolors )

    imsize = size( im(:,:,1) );
    displaycolors_hsv = rgb2hsv( displaycolors );
    
    imMultichannelOverlay = zeros( [imsize, 3] );
    
    cur_im_hsv = zeros( [imsize, 3] );
    for i = 1:size(im,3)
        curGrayImage = im(:,:,i);
        cur_im_hsv(:,:,1) = displaycolors_hsv(i,1);
        cur_im_hsv(:,:,2) = displaycolors_hsv(i,2);
        cur_im_hsv(:,:,3) = curGrayImage;  
        cur_im_rgb = hsv2rgb( cur_im_hsv );
        imMultichannelOverlay = imMultichannelOverlay + cur_im_rgb;
    end
    imMultichannelOverlay( imMultichannelOverlay > 1 ) = 1;
    
end