function [im_rgb] = blend_multichannel_image( im, displaycolors, blend_mode )
    
    if ~exist( 'blend_mode', 'var' )
        blend_mode = 'add';
    end
    
    if strcmpi( blend_mode , 'add' )
        
        % add their rgb images
        displaycolors_hsv = rgb2hsv( displaycolors );
        im_rgb = zeros( [size(im,1), size(im,2), 3] );
        cur_im_hsv = zeros( [size(im,1), size(im,2), 3] );
        for i = 1:size(im,3)
            curGrayImage = im(:,:,i);
            cur_im_hsv(:,:,1) = displaycolors_hsv(i,1);
            cur_im_hsv(:,:,2) = displaycolors_hsv(i,2);
            cur_im_hsv(:,:,3) = curGrayImage;  
            cur_im_rgb = hsv2rgb( cur_im_hsv );
            im_rgb = im_rgb + cur_im_rgb;
            im_rgb( im_rgb > 1 ) = 1;
        end
        im_rgb = im2uint8( im_rgb );
        
    elseif strcmpi( blend_mode , 'average' )
        
        % display average of the RGB images of all the channels
        displaycolors_hsv = rgb2hsv( displaycolors );
        im_rgb = zeros( [size(im,1), size(im,2), 3] );
        cur_im_hsv = zeros( [size(im,1), size(im,2), 3] );
        for i = 1:size(im,3)
            curGrayImage = im(:,:,i);
            cur_im_hsv(:,:,1) = displaycolors_hsv(i,1);
            cur_im_hsv(:,:,2) = displaycolors_hsv(i,2);
            cur_im_hsv(:,:,3) = curGrayImage;  
            cur_im_rgb = hsv2rgb( cur_im_hsv );
            im_rgb = im_rgb + cur_im_rgb / size(im,3);
        end
        im_rgb = im2uint8( im_rgb );
        
    else
        
        error( 'ERROR: Invalide blend_mode' );

    end
    
end

