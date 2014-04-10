function [ imCorrected ] = CorrectConfocalDepthSignalVariation_Opening( im, maxObjDiameter, flagDebugMode )

    if ~exist('flagDebugMode', 'var')
        flagDebugMode = false;
    end

    % compute means in each slice
    slicemean = zeros(1, size(im,3));
    slicestd = zeros(1, size(im,3));

    imBgndEstimate = imopen(im, ones(round(maxObjDiameter) * [1, 1]) );
    for z = 1:size( im, 3 )       
       
        imSlice = imBgndEstimate(:,:,z);
        slicemean(z) = mean( imSlice(:) );
        slicestd(z) = std( imSlice(:) );
        
    end

    if flagDebugMode
        imseriesshow(imBgndEstimate);
        figure, plot( slicemean );
    end
    
    % pick reference/standard mean
    [refmean, refSliceInd] = max( slicemean );
    refstd = slicestd(refSliceInd); 
    
    % corrrect each slice so that it has the mean equal to reference-mean
    % that we just picked
    imCorrected = im;
    
    for z = 1:size( im, 3 )       
        if slicestd(z) > eps
            imCorrected(:,:,z) = refmean + refstd * ((im(:,:,z) - slicemean(z)) / slicestd(z));
        end
    end
    
end