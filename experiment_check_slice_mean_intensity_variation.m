function experiment_check_slice_mean_intensity_variation(imageData, metadata)
%% check if slices are getting darker as we go deeper in z

    % normalize
%     for chid = 1:metadata.numChannels        
%         imageData{chid} = mat2gray(imageData{chid}) * 4096;
%     end
    
    % estimate foreground mask
    sbrCutoff = 2.0;
    diskKernelRadius = 30;
    resizeFactor = 1/2;
    
    krnlMax = streldisknd( round(diskKernelRadius * resizeFactor ./ metadata.voxelSpacing(1:2)) );
    imForegroundMask = cell(size(imageData));
    for chid = 1:metadata.numChannels
    
        imAdjusted = matitk( 'FMEDIAN', [1, 1, 1], imageData{chid} );
        
        imAdjusted = imresizend(imAdjusted, resizeFactor);
        imLocalBackground = imopen(imAdjusted, krnlMax);
        imSignalToBackgroundRatio = imAdjusted ./ (eps + imLocalBackground);
        imForegroundMask{chid} = imSignalToBackgroundRatio > sbrCutoff;
        imForegroundMask{chid} = imresizend(imForegroundMask{chid}, 1/resizeFactor, ...
                                            'interpolationMethod', 'nearest');
        
        imseriesmaskshow(imageData{chid}, imForegroundMask{chid});
        
    end
    
    % compute intensity stats accross slices
    sliceStats = cell(1, metadata.numChannels);
    for chid = 1:metadata.numChannels        
        for z = 1:metadata.volSize(3) 
            
            imSlice = imageData{chid}(:,:,z);   
            imMask = imForegroundMask{chid}(:,:,z);
            
            stats.max_cell = max(imSlice(imMask > 0));
            stats.mean_cell = mean(imSlice(imMask > 0));
            stats.std_cell = std(imSlice(imMask > 0));            

            stats.max_slice = max(imSlice(:));
            stats.mean_slice = mean(imSlice(:));
            stats.std_slice = std(imSlice(:));            
            
            sliceStats{chid}(z) = stats;
            
       end       
    end                  

   % generate plot
   numStd = 1;
   chColors = [0, 1, 0; 1 0 0];

   figure;
   hold on;
   zvals = ((1:metadata.volSize(3))-0.5) * metadata.voxelSpacing(3);
   for chid = 1:metadata.numChannels
       plot(zvals, [sliceStats{chid}.max_cell], 'Color', chColors(chid,:));
       ylim([0, 4096]);
       grid on;
   end
   hold off;   
   suptitle('Max cell intensity');
   
   figure;
   hold on;
   zvals = ((1:metadata.volSize(3))-0.5) * metadata.voxelSpacing(3);
   for chid = 1:metadata.numChannels
       mu = [sliceStats{chid}.mean_cell];
       lb = mu - numStd * [sliceStats{chid}.std_cell];
       ub = mu + numStd * [sliceStats{chid}.std_cell];
       
       subplot(1, metadata.numChannels, chid);
       signalHighLowPlot(zvals, mu, lb, ub, 'lineColor', chColors(chid,:));
       ylim([0, 4096]);
       grid on;
   end
   hold off;   
   suptitle('Mean cell intensity');
   

   figure;
   hold on;
   zvals = ((1:metadata.volSize(3))-0.5) * metadata.voxelSpacing(3);
   for chid = 1:metadata.numChannels
       mu = [sliceStats{chid}.mean_slice];
       lb = mu - numStd * [sliceStats{chid}.std_slice];
       ub = mu + numStd * [sliceStats{chid}.std_slice];
       
       subplot(1, metadata.numChannels, chid);
       signalHighLowPlot(zvals, mu, lb, ub, 'lineColor', chColors(chid,:));
       ylim([0, 4096]);
       grid on;
   end
   hold off;   
   suptitle('Mean slice intensity');
   
end

