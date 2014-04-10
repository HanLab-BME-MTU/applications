function experiment_check_slice_mean_intensity_variation( imageData, metadata )
%% check if slices are getting darker as we go deeper in z

    slice_mean_stats = zeros( metadata.numTimePoints, metadata.volSize(3) );
    slice_std_stats = zeros( metadata.numTimePoints, metadata.volSize(3) );
    slice_max_stats = zeros( metadata.numTimePoints, metadata.volSize(3) );
    figure;
    hold all;
    subplot(3,1,1), hold all; title( 'Mean image intensity in each slice' );
    subplot(3,1,2), hold all; title( 'Standard deviation image intensity in each slice' );
    subplot(3,1,3), hold all; title( 'Maximum image intensity in each slice' );
    for t = 1:metadata.numTimePoints      
       for z = 1:metadata.volSize(3) 
           imslice = imageData{t,1}(:,:,z);
           slice_mean_stats(t,z) = mean( imslice(:) );
           slice_std_stats(t,z) = std( imslice(:) );
           slice_max_stats(t,z) = max( imslice(:) );
       end                  
       subplot(3,1,1), plot( 1:metadata.volSize(3), slice_mean_stats(t,:) );
       subplot(3,1,2), plot( 1:metadata.volSize(3), slice_std_stats(t,:) );
       subplot(3,1,3), plot( 1:metadata.volSize(3), slice_max_stats(t,:) );
    end
    hold off;
    
end