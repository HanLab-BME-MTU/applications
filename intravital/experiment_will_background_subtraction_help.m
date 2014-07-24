function experiment_will_background_subtraction_help( imageData, metadata )
%% check if a background subtraction algorithm may help
%
% Essentially randomly shoot a bunch of rays through the 3D volume in 
% the time axis and visualize them to see if background subtraction can help

    numPoints = 5;    
    randPoints = floor( rand(numPoints,3) * diag( metadata.volSize ) );    
    imageData4DMat = cat( 4, imageData{:,1} ); 
    
    figure;
    
    for i = 1:numPoints
        subplot( numPoints, 1, i );        
            plot( 1:metadata.numTimePoints, squeeze( imageData4DMat(randPoints(i,1),randPoints(i,2),randPoints(i,3),:) ), 'r+-', 'LineWidth', 2.0 );
    end

end