function [cellColocStats] = ComputeDNADamageColocalizationMeasures(im53BP1, imColoc, spacing, ...
                                                                   imLabelCellSeg, cellStats, ...
                                                                   imColocSeg, varargin)
                                            
    p = inputParser;

    p.addParamValue('maxObjectRadius', 20, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('kernelDimensions', 2, @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x <= ndims(imDrug)));

    p.addParamValue('maxNeighDist', 50, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('numNeighDistLevels', 5, @(x) (isscalar(x) && isnumeric(x)));
    
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse(varargin{:});
          
    PARAMETERS = p.Results;
    
    numCells = numel(cellStats);
    volVoxel = prod(spacing);
    imdims = ndims(im53BP1);
    
    % pre-compute stuff used to compute measurements
    ptPixelInd = (1:numel(imColoc))';
    ptPixelLocPhysp = ind2submat(size(imColoc), ptPixelInd) * diag(spacing);
    
        % compute background corrected macrophage image
        imColocAdjusted = imColoc; %mat2gray(imColocAdjusted); % how well calibrated are the intensities?
        imColocAdjusted = matitk('FMEDIAN', round(2 * min(spacing) ./ spacing), double(imColocAdjusted));
        krnlMax = ones( 2 * round(PARAMETERS.maxObjectRadius ./ spacing(1:PARAMETERS.kernelDimensions)) + 1 ); % flat kernel is much more faster
        imColocAdjusted = imtophat(imColocAdjusted, krnlMax);
        
        % compute Coloc distance map
        imColocDistMap = bwdistsc(imColocSeg, spacing);

        % labeled cell boundary mask
        imLabelCellBoundarySeg = bwperim(imLabelCellSeg > 0) .* imLabelCellSeg;

    % compute measurmets for each cell    
    fprintf( '\nProgress: \n' );
    
    last_percent_done = 0;
    numPrint = 0;
    
    for cellId = 1:numCells
        
        curCellPixelIdxList = cellStats(cellId).PixelIdxList;

        ptCurCellPixelLoc = ind2submat(size(imLabelCellSeg), curCellPixelIdxList);
        ptCellMinCoord = min(ptCurCellPixelLoc);
        ptCellMaxCoord = max(ptCurCellPixelLoc);
        
        cellColocStats(cellId).cellId = cellId;
        
        % distance to nearest pixel in coloc segmentation mask
        cellColocStats(cellId).distanceToSegMask = min( imColocDistMap(curCellPixelIdxList) );

        % amount of macrophage inside cell nucleus
        cellColocStats(cellId).segPixelCount.InsideNucleus = sum(imColocSeg(curCellPixelIdxList) > 0);
        cellColocStats(cellId).segPixelVol.InsideNucleus = sum(imColocSeg(curCellPixelIdxList) > 0) * volVoxel;
        cellColocStats(cellId).segPixelFraction.InsideNucleus = mean(imColocSeg(curCellPixelIdxList) > 0);

        % macrophage intensity stats
        cellColocStats(cellId).intensityStats.insideNucleus = computeIntensityStats(imColocAdjusted, curCellPixelIdxList);

        curCellBBoxCropInd = cell(1,imdims);
        for i = 1:imdims
           minVal = max(1, ptCellMinCoord(i) - PARAMETERS.maxNeighDist);
           maxVal = min(size(imColoc,i), ptCellMaxCoord(i) + PARAMETERS.maxNeighDist);
           curCellBBoxCropInd{i} = minVal:maxVal;
        end

        imCurCellColocCropped = imColocAdjusted(curCellBBoxCropInd{:});
        imCurCellColocSegCropped = imColocSeg(curCellBBoxCropInd{:});
        imCurCellMaskCropped = (imLabelCellSeg(curCellBBoxCropInd{:}) == cellId);
        imCurCellDistMap = bwdistsc(imCurCellMaskCropped, spacing);

        binEdges = linspace(0, PARAMETERS.maxNeighDist, PARAMETERS.numMacrophageNeighDistLevels+1);

        for bid = 1:numel(binEdges)-1

            flagCurBin = (imCurCellDistMap > 0 & imCurCellDistMap <= binEdges(bid+1));
            ptCurBinPixelInd = find(flagCurBin);

            % stats from segmented coloc mask
            curSegPixelCount = sum( imCurCellColocSegCropped(ptCurBinPixelInd) > 0 );
            curSegPixelFraction = mean( imCurCellColocSegCropped(ptCurBinPixelInd) > 0 );

            strField = sprintf('segPixelCount.outsideNucleus.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount;', strField));

            strField = sprintf('segPixelCount.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount + cellColocStats(cellId).segPixelCount.InsideNucleus;', strField));

            strField = sprintf('segPixelVol.outsideNucleus.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount * volVoxel;', strField));

            strField = sprintf('segPixelVol.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = (curSegPixelCount * volVoxel) + cellColocStats(cellId).segPixelVol.InsideNucleus;', strField));

            strField = sprintf('segPixelFraction.outsideNucleus.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = curSegPixelFraction;', strField));

            strField = sprintf('segPixelFraction.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = mean( imCurCellColocSegCropped( imCurCellDistMap <= binEdges(bid+1) ) > 0 );', strField));
            
            % intensity stats within ribbons
            strField = sprintf('intensityStats.outsideNucleus.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellColocCropped, ptCurBinPixelInd);', strField));

            strField = sprintf('intensityStats.within_%d_um', binEdges(bid+1));
            eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellColocCropped, find(imCurCellDistMap <= binEdges(bid+1)));', strField));

        end
            
        % update progress/status
        percent_done = round(100*cellId/numCells);       
        
        if percent_done > last_percent_done
            fprintf( '%.2d%%  ', percent_done );
            last_percent_done = percent_done;
            numPrint = numPrint + 1;
            if mod( numPrint, 10 ) == 0
               fprintf( '\n' ); 
            end
        end        
        
    end
    
end
              
function istats = computeIntensityStats(imInput, pixelIdxList)

    pixelIntensities = imInput(pixelIdxList);
    
    istats.total = sum( pixelIntensities );

    istats.mean = mean( pixelIntensities );
    istats.std = std( pixelIntensities );
    istats.skewness = skewness( pixelIntensities );
    istats.kurtosis = kurtosis( pixelIntensities );

    istats.median = median( pixelIntensities );
    istats.mad = mad( pixelIntensities );
    istats.iqr = iqr( pixelIntensities );
    
end
                                                