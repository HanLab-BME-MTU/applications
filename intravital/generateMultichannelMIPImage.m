function [imResult] = generateMultichannelMIPImage(imageData, channelColorMap, spacing, gapSize)

    if ~exist( 'spacing', 'var' )
        spacing = ones(1,3);
    end

    if ~exist( 'gapSize', 'var' )
        gapSize = 10;
    end

    imStackMIP = cell(1,3);
    
    for dim = 1:3

        imStackMIP{dim} = [];

        for chid = 1:numel(imageData)

            imCurMIP = mat2gray(squeeze(max(imageData{chid}, [], dim)));

            if dim < 3
                outRows = size(imCurMIP, 1);
                outCols = round( size(imCurMIP, 2) * spacing(3) / spacing(dim) );
                imCurMIP = imresize(imCurMIP, [outRows, outCols]);
            end

            if dim == 1
                imCurMIP = imCurMIP';
            end

            imStackMIP{dim} = cat(3, imStackMIP{dim}, imCurMIP);

        end

    end
    
    
    if isempty( channelColorMap )

        imResult = [];
        
        for chid = 1:numel(imageData)    

            imCurResult = imStackMIP{3}(:,:,chid);
            imCurResult(end+1:end+gapSize, :) = 1;
            imCurResult = cat(1, imCurResult, imStackMIP{1}(:,:,chid));
            imCurResult(:,end+1:end+gapSize) = 1;
            imCurResult = cat(2, imCurResult, padarray(imStackMIP{2}(:,:,chid), [size(imCurResult,1)-size(imStackMIP{2},1), 0], 1.0, 'post'));

            if chid < numel(imageData)
                imCurResult(:,end+1:end+(2*gapSize)) = 1;
            end
            
            imResult = cat(2, imResult, imCurResult);

        end

        imResult = padarray(imResult, [gapSize, gapSize, 0], 1.0);
        
    else
        
        imResult = imStackMIP{3};
        imResult(end+1:end+gapSize,:,:) = 1;
        imResult = cat(1, imResult, imStackMIP{1});
        imResult(:,end+1:end+gapSize,:) = 1;
        imResult = cat(2, imResult, padarray(imStackMIP{2}, [size(imResult,1)-size(imStackMIP{2},1), 0, 0], 1.0, 'post'));

        imResult = genMultiChannelOverlay(imResult, channelColors);
        imResult = padarray(imResult, [gapSize, gapSize, 0], 1.0);
        
    end    
    
end