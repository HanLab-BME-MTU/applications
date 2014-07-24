function WriteCellSnapshotImages(imageData, imCellMask, cellStats, cellId, imageOutputDir)

    curCellCentroid = cellStats.Centroid;
    curCellBoundingBox = cellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), 70] );
    szOutputImage = [100, 100];

    % crop images
    subinds = cell(1,3);
    imsize = size(imageData{1});
    for i = 1:2

        xi = round(curCellCentroid(3-i) - 0.5 * curCellDisplaySize);

        xi_low = xi;
        if xi_low < 1 
            xi_low = 1;
        end

        xi_high = xi + curCellDisplaySize - 1;
        if xi_high > imsize(i)
            xi_high = imsize(i);
        end

        subinds{i} = xi_low:xi_high;

    end     
    subinds{3} = round(curCellCentroid(3));

    imCurCellMidSliceAllChannel = [];
    for j = 1:3
        channelDisplayRange = ComputeImageDynamicRange( imageData{j}, 98.0 );
        imCurChannelCropped = mat2gray( imageData{j}( subinds{:} ), channelDisplayRange );
        imCurChannelCropped = imresize( imCurChannelCropped, szOutputImage );
        imCurCellMidSliceAllChannel = cat(3, imCurCellMidSliceAllChannel, imCurChannelCropped );
    end

    imCurCellSegBndMidSliceCropped = imresize( bwperim( imCellMask( subinds{:} ) ), szOutputImage, 'nearest' );

    imCurCellMIP = mat2gray( max( imageData{1}(subinds{1:2}, :) .* imCellMask(subinds{1:2}, :), [], 3) );
    imCurCellMIP = imresize(imCurCellMIP, szOutputImage);

    % write images
    imwrite( genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   

    channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
    imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   

    imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceAllChannel_%.3d.png', cellId)), 'png' );   

    imwrite(imCurCellMIP, fullfile(imageOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   

end