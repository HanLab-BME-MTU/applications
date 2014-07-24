function [ imCorrected ] = CorrectConfocalDepthSignalVariation( im, maxObjDiameter, flagDebugMode )

    if ~exist('flagDebugMode', 'var')
        flagDebugMode = false;
    end
    
    % compute means in each slice
    slicemean = zeros( 1, size(im,3) );
    imSmooth = zeros( size(im) );
    
    smoothingKernel = fspecial('average', 2*maxObjDiameter);
    for z = 1:size( im, 3 )       
        
        imSlice = imfilter( im(:,:,z), smoothingKernel );
        slicemean(z) = max( imSlice(:) );        
        imSmooth(:,:,z) = imSlice;
        
    end

    if flagDebugMode
        imseriesshow(imSmooth);
        figure, plot( slicemean );
    end
    
    % pick reference/standard mean
    refmean = max( slicemean );
    
    % corrrect each slice so that it has the mean equal to reference-mean
    % that we just picked
    imCorrected = im;
    for z = 1:size( im, 3 )       
        %imCorrected(:,:,z) = im(:,:,z) * (refmean / (slicemean(z) + eps));
        %imCorrected(:,:,z) = im(:,:,z) + (refmean - slicemean(z));
        imCorrected(:,:,z) = mat2gray(im(:,:,z)) * refmean;
    end
    
end