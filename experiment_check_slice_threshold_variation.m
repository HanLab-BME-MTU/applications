function experiment_check_slice_threshold_variation( imageData, imageDataMIP, metadata )
%% compute minimum error thresholds for each slice and see how they vary    
    th_level1 = zeros( metadata.numTimePoints, metadata.volSize(3) );    
    th_level2 = zeros( metadata.numTimePoints, metadata.volSize(3) );    
    imageDataMask_level1 = imageData;
    imageDataMask_level2 = imageData;
    figure; hold all;
    for t = 1:metadata.numTimePoints
        for z = 1:metadata.volSize(3)
            imslice = imageData{t,1}(:,:,z);
            th_level1(t,z) = th_minerror( imslice, 2000 );        
            imageDataMask_level1{t,1}(:,:,z) = imslice >= th_level1(t,z);

            th_level2(t,z) = th_minerror( imslice( imslice >= th_level1(t,z) ), 2000 );        
            imageDataMask_level2{t,1}(:,:,z) = imslice >= th_level2(t,z);
        end
        plot( 1:metadata.volSize(3), th_level1(t,:) );
    end

    th_mip_level1 = zeros( metadata.numTimePoints, metadata.numChannels );
    th_mip_level2 = zeros( metadata.numTimePoints, metadata.numChannels );
    imageDataMIPMask_level1 = imageDataMIP;
    imageDataMIPMask_level2 = imageDataMIP;    
    for c = 1:metadata.numChannels 
        for t = 1:metadata.numTimePoints
            imslice = imageDataMIP{c}(:,:,t);

            th_mip_level1(t) = th_minerror( imslice, 2000 );        
            imageDataMIPMask_level1{t,1}(:,:,z) = imslice >= th_level1(t);

            th_mip_level2(t) = th_minerror( imslice( imslice >= th_level1(t) ), 2000 );        
            imageDataMIPMask_level2{t,1}(:,:,z) = imslice >= th_level2(t,z);            
        end
    end 
end