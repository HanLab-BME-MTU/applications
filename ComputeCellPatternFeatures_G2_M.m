function [classNameList, featureStruct, featureVec, featureNameList, flagLearnFeature] = ComputeCellFeatures_G2_M( imageData, imCellSegMask, cellStats, flagAnnotated )

    featureStruct = [];
    featureVec = [];
    featureNameList = [];
    flagLearnFeature = [];
    
    if ~exist( 'flagAnnotated', 'var' )
        flagAnnotated = true;
    end
    
    % assign class label
    classNameList = { 'G2', 'M' };
    
    classNameList = sort(classNameList);
    if nargin == 0
        return;
    end
    
    if flagAnnotated
        
        switch(cellStats.cellPatternType)

            case { 'Mono_G2', 'Multi_G2'}     

                featureStruct.label.className = 'G2';

            case { 'Mono_Prophase', 'Multi_Prophase' }     

                featureStruct.label.className = 'M';

            otherwise

                return;

        end    
    
        featureStruct.label.classid = strmatch( featureStruct.label.className, classNameList );
        
    else
        
        featureStruct.label.className = 'G2'; % dummy  label
        
    end
    
    % note down some cell identification info for debugging -- should be
    % ignored in classfication
    if flagAnnotated
        
        featureStruct.cellInfo.annotationFile = cellStats.annotationFile;
        featureStruct.cellInfo.cellId = cellStats.cellId;
        featureStruct.cellInfo.flagIsOnBorder = cellStats.flagIsOnBorder;
        featureStruct.cellInfo.annotatedCellPatternType = cellStats.cellPatternType;
        featureStruct.cellInfo.annotatedCellPatternId = cellStats.cellPatternId;    

    end
    
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

    imCellCropped = cell(size(imageData));
    for i = 1:numel(imageData)
        imCellCropped{i} = imageData{i}( subinds{1:3} );
    end
    
    imCellSegCropped = imCellSegMask( subinds{1:3} );
    
    % intensity statistics of cell pixels in blue and green channels
    imCellMaskEroded = imerode(imCellSegCropped, streldisknd([1,1,1]));
    for i = 1:1        
        
        cellPixelIntensities = imCellCropped{i}(imCellSegCropped > 0);                    
        
        % global intensity
        featureStruct.intensity.std(i) = std( cellPixelIntensities );
        featureStruct.intensity.skewness(i) = skewness( cellPixelIntensities );
        featureStruct.intensity.kurtosis(i) = kurtosis( cellPixelIntensities );
        featureStruct.intensity.entropy(i) = entropy( cellPixelIntensities );
        
        % local intensity        
        imLocalStd = stdfilt( imCellCropped{i} );
        featureStruct.intensity.LocalStddev = mean2( imLocalStd(imCellMaskEroded > 0) );
        
        imLocalEntropy = entropyfilt( imCellCropped{i} );
        featureStruct.intensity.LocalEntropy = mean2( imLocalEntropy(imCellMaskEroded > 0) );
        
    end
    
    % GLCM texture features
    disp = [2, 4];
    OffsetDictionary = [0  1 0;  1 0 0;  1  1 0; -1  1 0; ...
                        0  1 1;  1 0 1;  1  1 1; -1  1 1; ...
                        0 -1 1; -1 0 1; -1 -1 1;  1 -1 1; ...
                        0  0 1];
               
    Offsets = [];
    for i = 1:numel(disp)
        Offsets = [Offsets; OffsetDictionary * disp(i)];
    end
    Offsets(:,3) = Offsets(:,3) / 2.0;
    
    ImageIntensityRange = ComputeImageIntensityRange( imCellCropped{1} );
    
    [GLCMS] = graycomatrixnd( imCellCropped{1}, ...
                              'ROIMask', imCellCropped{1}, ...  
                              'Offset', Offsets, ...
                              'NumLevels', 64, ...
                              'GrayLimits', ImageIntensityRange, ...
                              'Symmetric', true );
     
    featureStruct.texture.glcm = graycopropsext( sum(GLCMS, 3) );
    
    % Generate concentric band features
    
    % generate feature vector
    [ featureVec , featureNameList ] = ConvertFeatureStructToFeatureVec( featureStruct );
    
    % flag on features that shouldbe used in machine learning
    flagLearnFeature = false(numel(featureNameList),1);
    featureExcludeNameList = {'label.classid', 'cellInfo'};
    for i = 1:numel(featureExcludeNameList)
        flagLearnFeature = flagLearnFeature | strncmp( featureNameList, featureExcludeNameList{i}, numel(featureExcludeNameList{i}) );
    end
    flagLearnFeature = ~flagLearnFeature;
    
end


% TO BE IMPLEMENTED
% %% difference of intensity statistics of nuclear channel in three concentric bands
% cellBBox = { min(yind):max(yind), min(xind):max(xind), min(zind):max(zind) };
% imCellCropped = imageData{1}( cellBBox );
% imCellMaskCropped = imCellSegMask{1}( cellBBox );
    
