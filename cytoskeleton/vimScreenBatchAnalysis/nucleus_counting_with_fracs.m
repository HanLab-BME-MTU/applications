function  nucleus_count = ...
    nucleus_counting_with_fracs(OtsuRosin_Segment, ...
    currentImg,NucleusSegmentationChannelOutputDir, iChannel, iFrame,display_save_flag)

[label_img, nucleus_count] = bwlabel(OtsuRosin_Segment);

if(display_save_flag>0)
    h = figure;
    subplot(211);
    imagesc(currentImg);axis off;axis equal;colormap(gray);colorbar('EastOutside');freezeColors;
    subplot(212);    
    imagesc(label_img);colormap(jet); axis off;  axis equal;
    % permute colors for better display
    C = colormap;
    new_C = [1 1 1;C(randperm(size(C,1)-1),:)];
    colormap(new_C);colorbar('EastOutside');
    title(['Number of nucleus:', num2str(nucleus_count)]);
    saveas(h,[NucleusSegmentationChannelOutputDir,filesep, ...
        'nucleus_count_im_c', num2str(iChannel), '_f',num2str(iFrame),'.fig']);
    saveas(h,[NucleusSegmentationChannelOutputDir,filesep, ...
        'nucleus_count_im_c', num2str(iChannel), '_f',num2str(iFrame),'.tif']);
    close(h);
    
    h = figure;    
    imagesc(label_img);colormap(jet); axis off;  axis equal;
    % permute colors for better display
    C = colormap;
    new_C = [1 1 1;C(randperm(size(C,1)-1),:)];
    colormap(new_C);colorbar('EastOutside');
    title(['Number of nucleus:', num2str(nucleus_count)]);
    saveas(h,[NucleusSegmentationChannelOutputDir,filesep, ...
        'nucleus_count_c', num2str(iChannel), '_f',num2str(iFrame),'.fig']);
    saveas(h,[NucleusSegmentationChannelOutputDir,filesep, ...
        'nucleus_count_c', num2str(iChannel), '_f',num2str(iFrame),'.tif']);
end

      