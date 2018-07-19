function [featureStruct, featureVec, featureNameList, flagLearnFeature] = ComputeCellPatternFeatures_PreG1_G1_S_G2M( imageData, imCellSegMask, cellStats )

    flagDebug = true;
    flagApplyLog = true;
    
    featureStruct = [];
    featureVec = [];
    featureNameList = [];
    flagLearnFeature = [];
    
    % assign class label
    classNameList = { 'Pre-G1', 'G1', 'S', 'G2-M' };
    
    switch(cellStats.cellPatternType)
       
        case { 'Mono_Pre_G1', 'Multi_Pre_G1' } 

            featureStruct.label.className = 'Pre-G1';

        case { 'Mono_G1', 'Multi_G1' } 

            featureStruct.label.className = 'G1';

        case { 'Mono_S', 'Multi_S' } 

            featureStruct.label.className = 'S';

        case { 'Mono_G2', 'Multi_G2', 'Mono_Prophase', 'Multi_Prophase' }     

            featureStruct.label.className = 'G2-M';

        otherwise
        
            return;
            
    end             
        
    featureStruct.label.classid = strmatch( featureStruct.label.className, classNameList );

    % note down some cell identification info for debugging -- should be
    % ignored in classfication
    featureStruct.cellInfo.annotationFile = cellStats.annotationFile;
    featureStruct.cellInfo.cellId = cellStats.cellId;
    featureStruct.cellInfo.flagIsOnBorder = cellStats.flagIsOnBorder;
    featureStruct.cellInfo.annotatedCellPatternType = cellStats.cellPatternType;
    featureStruct.cellInfo.annotatedCellPatternId = cellStats.cellPatternId;    
    
    % crop a small 3D patch enclosing the cell
    cropBoxSize = 80;
    cellCentroid = cellStats.Centroid;
    cellBoundingBox = cellStats.BoundingBox;
    cellDisplaySize = max( [cellBoundingBox(4:5), cropBoxSize] );
    
    subinds = cell(1,3);
    imsize = size(imageData{1});
    for i = 1:2

        xi = round(cellCentroid(3-i) - 0.5 * cellDisplaySize);

        xi_low = xi;
        if xi_low < 1 
            xi_low = 1;
        end

        xi_high = xi + cellDisplaySize - 1;
        if xi_high > imsize(i)
            xi_high = imsize(i);
        end

        subinds{i} = xi_low:xi_high;

    end    
    subinds{3} = round(cellStats.BoundingBox(3):(cellStats.BoundingBox(3)+cellStats.BoundingBox(6)-1));

    imCellSegCropped = imCellSegMask( subinds{1:3} );

    % extract cell mask in each channel by local thresholding -- this sort
    % of makes it robust/invariant to any inaccurate shift correction
    % between channels
    imsize = size(imageData{1});
    
    imCellImageDataCropped = cell(1,3);
    imCellMaskCropped = cell(1,3);
    
    subIndsCellZStack = subinds;
    subIndsCellZStack{1} = 1:imsize(1);
    subIndsCellZStack{2} = 1:imsize(2);
    
    flagStateChannelVisibility = [ 0 0 ; 0 1 ; 1 1 ; 1 0 ];
    
    for i = 1:numel(imageData)        
        
        imCellImageDataCropped{i} = imageData{i}( subinds{:} );
        
        if i > 1
            
            imCellMaskCropped{i} = false( size(imCellImageDataCropped{i}) );
            
            curGlobalThresholdVal = thresholdOtsu( imageData{i}( subIndsCellZStack{:} ) );
            curLocalThreshVal = thresholdOtsu( imCellImageDataCropped{i} );
            
            flagOverlapTooLess = false;
            
            if curLocalThreshVal > 0.25 * curGlobalThresholdVal
                imCellMaskCropped{i} = imCellImageDataCropped{i} > curLocalThreshVal;
                imCellMaskCropped{i} = imCellMaskCropped{i} & imCellSegCropped > 0;
                regstats = regionprops( imCellMaskCropped{i}, { 'Area', 'PixelIdxList' } );
                if numel(regstats) >= 1
                    [maxArea, maxind] = max( [regstats.Area] );
                    imCellMaskCropped{i} = false( size(imCellImageDataCropped{i}) );
                    
                    overlapMagnitude = maxArea / numel(cellStats.PixelIdxList);
                    
                    % there must be a significant overlap                        
                    if overlapMagnitude > 0.10 
                        imCellMaskCropped{i}( regstats(maxind).PixelIdxList ) = true;
                    else
                        flagOverlapTooLess = true;                        
                    end
                    
                end
            end
            
            if flagDebug 
                if ~any( imCellMaskCropped{i}(:) ) 
                    flagChannelVisible = false;
                    strVisible = 'empty';
                else
                    flagChannelVisible = true;
                    strVisible = 'not empty';
                end
                
                if flagChannelVisible ~= flagStateChannelVisibility( featureStruct.label.classid, i - 1 ) 
                    
                    if ~flagChannelVisible && flagOverlapTooLess
                        fprintf( '\nWARNING: Magnitue of Overlap with Cell mask of channel %d is too less\n', i );                        
                    end
                    
                    fprintf( '\nWARNING: Cell mask of channel %d is %s for %s cell with id = %d\n', ...
                             i, strVisible, cellStats.cellPatternType, cellStats.cellId );
                end
            end
        else
            imCellMaskCropped{i} = imCellSegCropped;
        end
            
    end
    
    % median, mean, min, max, and std cell intensity in all channels
    featureStruct.intensity.median = zeros(1,numel(imageData));
    featureStruct.intensity.mean = zeros(1,numel(imageData));
    featureStruct.intensity.min = zeros(1,numel(imageData));
    featureStruct.intensity.max = zeros(1,numel(imageData));
    featureStruct.intensity.stddev = zeros(1,numel(imageData));
    
    for i = 1:numel(imageData)        
        curImageData = imCellImageDataCropped{i};
        curCellMask = imCellMaskCropped{i};
        
        if any( curCellMask(:) )            
            featureStruct.intensity.median(i) = median( curImageData(curCellMask) );
            featureStruct.intensity.mean(i) = mean( curImageData(curCellMask) );
            featureStruct.intensity.stddev(i) = std( curImageData(curCellMask) );
            featureStruct.intensity.min(i) = min( curImageData(curCellMask) );
            featureStruct.intensity.max(i) = max( curImageData(curCellMask) );
        end

%         featureStruct.intensity.median(i) = median( imageData{i}(cellStats.PixelIdxList) );
%         featureStruct.intensity.mean(i) = mean( imageData{i}(cellStats.PixelIdxList) );
%         featureStruct.intensity.stddev(i) = std( imageData{i}(cellStats.PixelIdxList) );
%         featureStruct.intensity.min(i) = min( imageData{i}(cellStats.PixelIdxList) );
%         featureStruct.intensity.max(i) = max( imageData{i}(cellStats.PixelIdxList) );
    end
    
    % compute meta intensity features
    infVal = 4096;
    chcombs = combnk(1:numel(imageData),2);
    metaFeatureNameList = { 'median', 'mean', 'min', 'max' };
    for i = 1:numel(metaFeatureNameList)       
        curFeatureArray = getfield( featureStruct.intensity, metaFeatureNameList{i} );
        for j = 1:size(chcombs,1)        
            curComb = chcombs(j,:);        
            curFeatureName = sprintf( '%sRatio_ch%d_by_%d', ...
                                      metaFeatureNameList{i}, ...
                                      curComb(2), curComb(1) );         
            curFeatureVal = curFeatureArray( curComb(2) ) / curFeatureArray( curComb(1) );    
            curFeatureVal = FixInvalidFeatureVal( curFeatureVal, infVal );
            featureStruct.intensity = setfield( featureStruct.intensity, ...
                                                curFeatureName, ...
                                                curFeatureVal );
        end
    end
    
    % generate feature vector
    [ featureVec , featureNameList ] = ConvertFeatureStructToFeatureVec( featureStruct );
    
    % flag on features that shouldbe used in machine learning
    flagLearnFeature = false(numel(featureNameList),1);
    featureExcludeNameList = {'label', 'cellInfo'};
    for i = 1:numel(featureExcludeNameList)
        flagLearnFeature = flagLearnFeature | strncmp( featureNameList, featureExcludeNameList{i}, numel(featureExcludeNameList{i}) );
    end
    flagLearnFeature = ~flagLearnFeature;
    
end

function [ featureValFixed ] = FixInvalidFeatureVal( featureVal, infVal )

    if isnan( featureVal )
        featureValFixed = 0;
    elseif isinf( featureVal )
        featureValFixed = infVal;
    else
        featureValFixed = featureVal;
    end
    
end