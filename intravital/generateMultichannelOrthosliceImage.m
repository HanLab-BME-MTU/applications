function [imResult] = generateMultichannelOrthosliceImage(imageData, pixelPos, channelColorMap, spacing, gapSize)

    if ~exist( 'spacing', 'var' )
        spacing = ones(1,3);
    end

    if ~exist( 'gapSize', 'var' )
        gapSize = 10;
    end
    
    imOrthoSlices = cell(1,3);

    for dim = 1:3

        imOrthoSlices{dim} = [];

        indImageSlicer = { ':', ':', ':' };
        indImageSlicer{dim} = pixelPos(dim);

        for chid = 1:numel(imageData)

            imCurOrthoSlice = mat2gray(squeeze(imageData{chid}(indImageSlicer{:})));

            if dim < 3
                outRows = size(imCurOrthoSlice, 1);
                outCols = round( size(imCurOrthoSlice, 2) * spacing(3) / spacing(dim) );
                imCurOrthoSlice = imresize(imCurOrthoSlice, [outRows, outCols]);
            end

            if dim == 1
                imCurOrthoSlice = imCurOrthoSlice';
            end

            imOrthoSlices{dim} = cat(3, imOrthoSlices{dim}, imCurOrthoSlice);

        end

    end
    
    imResult = imOrthoSlices{3};
    imResult(end+1:end+gapSize,:,:) = 1;
    imResult = cat(1, imResult, imOrthoSlices{1});
    imResult(:,end+1:end+gapSize,:) = 1;
    imResult = cat(2, imResult, padarray(imOrthoSlices{2}, [size(imResult,1)-size(imOrthoSlices{2},1), 0, 0], 1.0, 'post'));

    imResult = genMultiChannelOverlay(imResult, channelColorMap);
    imResult = padarray(imResult, [gapSize, gapSize, 0], 1.0);
    
end