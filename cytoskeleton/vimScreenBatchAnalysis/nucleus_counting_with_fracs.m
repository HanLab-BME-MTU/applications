function  nucleus_count = ...
    nucleus_counting_with_fracs(OtsuRosin_Segment, ...
    currentImg,NucleusSegmentationChannelOutputDir, iChannel, iFrame,display_save_flag)

[label_img, nucleus_count] = bwlabel(OtsuRosin_Segment);

if(display_save_flag>0)
    h = figure;
    subplot(121);
    imagesc(currentImg);colormap(gray);axis off;
    subplot(122);    
    imagesc(label_img);colorbar('EastOutside'); axis off;  
    % permute colors for better display
    C = colormap;
    new_C = C(randperm(size(C,1)),:);
    colormap(new_C);
    title(['Number of nucleus:', num2str(nucleus_count)]);
    saveas(h,[NucleusSegmentationChannelOutputDir,filesep, ...
        'nucleus_count_c', num2str(iChannel), '_f',num2str(iFrame),'tif']);
end

      