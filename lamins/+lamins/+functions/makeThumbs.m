function [ output_args ] = makeThumbs( MD, order )
%makeThumbs Make thumbnails
    if(nargin < 2)
        order = 1:4;
    end
    I = CellReader(MD.getReader());
    thumbs = vertcat(imadjust(horzcat(I{order(1),:})),imadjust(horzcat(I{order(2),:})),imadjust(horzcat(I{order(3),:})),imadjust(horzcat(I{order(4),:})));
    
    name = MD.getFilename();
    name = strrep(name,' ','_');
    name = strrep(name,'.mat','');

    %imwrite(im2uint8(thumbs),jet(256),'thumbs2.png');
    thumbs = imresize(thumbs,0.25);
    imwrite(im2uint8(thumbs),jet(256),['small_thumbs_' name '.png']);
    thumbs = imresize(thumbs,0.25);
    imwrite(im2uint8(thumbs),jet(256),['tiny_thumbs_' name ' .png']);
    thumbs = imresize(thumbs,0.25);
    imwrite(im2uint8(thumbs),jet(256),['mini_thumbs_' name '.png']);
    thumbs = imresize(thumbs,0.25);
    imwrite(im2uint8(thumbs),jet(256),['micro_thumbs_' name ' .png']);
end

