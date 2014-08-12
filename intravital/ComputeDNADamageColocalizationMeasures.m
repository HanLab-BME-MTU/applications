function [cellColocStats] = ComputeDNADamageColocalizationMeasures(imDrug, im53BP1, imMacrophage, spacing, ...
                                                                   imLabelCellSeg, cellStats, ...
                                                                   imDrugSeg, imMacrophageSeg, varargin)
                                            
    p = inputParser;

    p.addParamValue('kernelDimensions', 2, @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x <= ndims(imDrug)));
    
    p.addParamValue('maxMacrophageObjectRadius', 20, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('maxMacrophageNeighDist', 50, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('numMacrophageNeighDistLevels', 5, @(x) (isscalar(x) && isnumeric(x)));
    
    p.addParamValue('maxDrugObjectRadius', 20, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('maxDrugNeighDist', 50, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('numDrugNeighDistLevels', 5, @(x) (isscalar(x) && isnumeric(x)));

    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse(varargin{:});
          
    PARAMETERS = p.Results;
    
    numCells = numel(cellStats);
    volVoxel = prod(spacing);
    imdims = ndims(im53BP1);
    
    % pre-compute stuff used to compute measurements
    ptPixelInd = (1:numel(imMacrophage))';
    ptPixelLocPhysp = ind2submat(size(imMacrophage), ptPixelInd) * diag(spacing);
    
        % compute background corrected drug image
        imDrugAdjusted = imDrug; % mat2gray(imDrug); % how well calibrated are the intensities?
        imDrugAdjusted = matitk('FMEDIAN', round(2 * min(spacing) ./ spacing), double(imDrugAdjusted));
        krnlMax = ones( 2 * round(PARAMETERS.maxDrugObjectRadius ./ spacing(1:PARAMETERS.kernelDimensions)) + 1 );
        imDrugAdjusted = imtophat(imDrugAdjusted, krnlMax);
        
        % compute background corrected macrophage image
        imMacrophageAdjusted = imMacrophage; %mat2gray(imMacrophageAdjusted); % how well calibrated are the intensities?
        imMacrophageAdjusted = matitk('FMEDIAN', round(2 * min(spacing) ./ spacing), double(imMacrophageAdjusted));
        krnlMax = ones( 2 * round(PARAMETERS.maxMacrophageObjectRadius ./ spacing(1:PARAMETERS.kernelDimensions)) + 1 ); % flat kernel is much more faster
        imMacrophageAdjusted = imtophat(imMacrophageAdjusted, krnlMax);
        
        % compute macrophage distance map
        imMacrophageDistMap = bwdistsc(imMacrophageSeg, spacing);

        % get drug pixel locations
        drugPixelInd = find(imDrugSeg > 0);
        ptDrugPixelLocPhysp = ind2submat(size(imDrugSeg), drugPixelInd) * diag(spacing);

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
        
        % measurements from macrophage channel
        
            % distance to nearest macrophage
            cellColocStats(cellId).macrophage.distanceToMacrophage = min( imMacrophageDistMap(curCellPixelIdxList) );

            % amount of macrophage inside cell nucleus
            cellColocStats(cellId).macrophage.segPixelCountInsideNucleus = sum(imMacrophageSeg(curCellPixelIdxList) > 0);
            cellColocStats(cellId).macrophage.segPixelVolInsideNucleus = sum(imMacrophageSeg(curCellPixelIdxList) > 0) * volVoxel;
            cellColocStats(cellId).macrophage.segPixelFractionInsideNucleus = mean(imMacrophageSeg(curCellPixelIdxList) > 0);
        
            % macrophage intensity stats
            cellColocStats(cellId).macrophage.intensityStats.insideNucleus = computeIntensityStats(imMacrophageAdjusted, curCellPixelIdxList);

            curCellBBoxCropInd = cell(1,imdims);
            for i = 1:imdims
               minVal = max(1, ptCellMinCoord(i) - PARAMETERS.maxMacrophageNeighDist);
               maxVal = min(size(imMacrophage,i), ptCellMaxCoord(i) + PARAMETERS.maxMacrophageNeighDist);
               curCellBBoxCropInd{i} = minVal:maxVal;
            end
            
            imCurCellMacrophageCropped = imMacrophageAdjusted(curCellBBoxCropInd{:});
            imCurCellMacrophageSegCropped = imMacrophageSeg(curCellBBoxCropInd{:});
            imCurCellMaskCropped = (imLabelCellSeg(curCellBBoxCropInd{:}) == cellId);
            imCurCellDistMap = bwdistsc(imCurCellMaskCropped, spacing);

            binEdges = linspace(0, PARAMETERS.maxMacrophageNeighDist, PARAMETERS.numMacrophageNeighDistLevels+1);
            
            for bid = 1:numel(binEdges)-1

                flagCurBin = (imCurCellDistMap > 0 & imCurCellDistMap <= binEdges(bid+1));
                ptCurBinPixelInd = find(flagCurBin);
                
                % stats from segmented macrophage mask
                curSegPixelCount = sum( imCurCellMacrophageSegCropped(ptCurBinPixelInd) > 0 );
                
                strField = sprintf('macrophage.segPixelCountOutsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount;', strField));

                strField = sprintf('macrophage.segPixelVolOutsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount * volVoxel;', strField));

                strField = sprintf('macrophage.segPixelCount.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount + cellColocStats(cellId).macrophage.segPixelCountInsideNucleus;', strField));
                
                strField = sprintf('macrophage.segPixelVolwithin_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = (curSegPixelCount * volVoxel) + cellColocStats(cellId).macrophage.segPixelVolInsideNucleus;', strField));

                % intensity stats within ribbons
                strField = sprintf('macrophage.intensityStats.outsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellMacrophageCropped, ptCurBinPixelInd);', strField));

                strField = sprintf('macrophage.intensityStats.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellMacrophageCropped, find(imCurCellDistMap <= binEdges(bid+1)));', strField));
                
            end
            
        % measurments from drug channel 

            % amount of drug inside cell nucleus
            cellColocStats(cellId).drug.segPixelCountInsideNucleus = sum(imDrugSeg(curCellPixelIdxList) > 0);
            cellColocStats(cellId).drug.segPixelVolInsideNucleus = sum(imDrugSeg(curCellPixelIdxList) > 0) * volVoxel;
            cellColocStats(cellId).drug.segPixelFractionInsideNucleus = mean(imDrugSeg(curCellPixelIdxList) > 0);

            % drug intensity stats
            cellColocStats(cellId).drug.intensityStats.insideNucleus = computeIntensityStats(imDrugAdjusted, curCellPixelIdxList);
            
            curCellBBoxCropInd = cell(1,imdims);
            for i = 1:imdims
               minVal = max(1, ptCellMinCoord(i) - PARAMETERS.maxDrugNeighDist);
               maxVal = min(size(imDrug,i), ptCellMaxCoord(i) + PARAMETERS.maxDrugNeighDist);
               curCellBBoxCropInd{i} = minVal:maxVal;
            end
            
            imCurCellDrugCropped = imDrugAdjusted(curCellBBoxCropInd{:});
            imCurCellDrugSegCropped = imDrugSeg(curCellBBoxCropInd{:});
            imCurCellMaskCropped = (imLabelCellSeg(curCellBBoxCropInd{:}) == cellId);
            imCurCellDistMap = bwdistsc(imCurCellMaskCropped, spacing);

            binEdges = linspace(0, PARAMETERS.maxDrugNeighDist, PARAMETERS.numDrugNeighDistLevels+1);
            
            for bid = 1:numel(binEdges)-1

                flagCurBin = (imCurCellDistMap > 0 & imCurCellDistMap <= binEdges(bid+1));
                ptCurBinPixelInd = find(flagCurBin);
                
                % stats from segmented macrophage mask
                curSegPixelCount = sum( imCurCellDrugSegCropped(ptCurBinPixelInd) > 0 );
                
                strField = sprintf('drug.segPixelCountOutsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount;', strField));

                strField = sprintf('drug.segPixelVolOutsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount * volVoxel;', strField));

                strField = sprintf('drug.segPixelCount.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount + cellColocStats(cellId).drug.segPixelCountInsideNucleus;', strField));
                
                strField = sprintf('drug.segPixelVolwithin_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = (curSegPixelCount * volVoxel) + cellColocStats(cellId).drug.segPixelVolInsideNucleus;', strField));

                % intensity stats within ribbons
                strField = sprintf('drug.intensityStats.outsideNucleus.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellDrugCropped, ptCurBinPixelInd);', strField));

                strField = sprintf('drug.intensityStats.within_%d_um', binEdges(bid+1));
                eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imCurCellDrugCropped, find(imCurCellDistMap <= binEdges(bid+1)));', strField));
                
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
    
%     % compute measurmets for each cell    
%     fprintf( '\nProgress: \n' );
%     
%     last_percent_done = 0;
%     numPrint = 0;
%     
%     for cellId = 1:numCells
%         
%         curCellPixelIdxList = cellStats(cellId).PixelIdxList;
% 
%         curCellBoundaryPixelIdxList = find(imLabelCellBoundarySeg == cellId);
%         ptCurCellBoundaryPixelLocPhysp = ind2submat(size(imLabelCellSeg), curCellBoundaryPixelIdxList) * diag(spacing);
%         ptCellMinCoord = min(ptCurCellBoundaryPixelLocPhysp);
%         ptCellMaxCoord = max(ptCurCellBoundaryPixelLocPhysp);
%         kdtree_ptcell = KDTreeSearcher(ptCurCellBoundaryPixelLocPhysp);          
%         
%         cellColocStats(cellId).cellId = cellId;
%         
%         % measurements from macrophage channel
%         
%             % distance to nearest macrophage
%             cellColocStats(cellId).macrophage.distanceToMacrophage = min( imMacrophageDistMap(curCellPixelIdxList) );
% 
%             % amount of macrophage inside cell nucleus
%             cellColocStats(cellId).macrophage.segPixelCountInsideNucleus = sum(imMacrophageSeg(curCellPixelIdxList) > 0);
%             cellColocStats(cellId).macrophage.segPixelVolInsideNucleus = sum(imMacrophageSeg(curCellPixelIdxList) > 0) * volVoxel;
%             cellColocStats(cellId).macrophage.segPixelFractionInsideNucleus = mean(imMacrophageSeg(curCellPixelIdxList) > 0);
%         
%             % macrophage intensity stats
%             cellColocStats(cellId).macrophage.intensityStats.insideNucleus = computeIntensityStats(imMacrophageAdjusted, curCellPixelIdxList);
%             
%             flagInValidNeighborhood = true(numel(imMacrophage), 1);
%             flagInValidNeighborhood(curCellPixelIdxList) = false;
%             for j = 1:ndims(imMacrophage)
%                 flagInValidNeighborhood(ptPixelLocPhysp(:,j) <= ptCellMinCoord(j) - PARAMETERS.maxMacrophageNeighDist | ...
%                                         ptPixelLocPhysp(:,j) >= ptCellMaxCoord(j) + PARAMETERS.maxMacrophageNeighDist) = false;
%             end
%             ptCurPixelLocPhysp = ptPixelLocPhysp(flagInValidNeighborhood, :);
%             ptCurPixelInd = ptPixelInd(flagInValidNeighborhood);
%             
%             [ptIdNearest, dist] = knnsearch(kdtree_ptcell, ptCurPixelLocPhysp);
%             flagValidDist = (dist > 0 & dist < PARAMETERS.maxMacrophageNeighDist);
%             ptCurPixelInd = ptCurPixelInd(flagValidDist);
%             dist = dist(flagValidDist);
%             
%             binEdges = linspace(0, PARAMETERS.maxMacrophageNeighDist, PARAMETERS.numMacrophageNeighDistLevels+1);
%             
%             for bid = 1:numel(binEdges)-1
% 
%                 flagCurBin = (dist >= binEdges(bid) & dist < binEdges(bid+1));
%                 ptCurBinPixelInd = ptCurPixelInd(flagCurBin);
%                 
%                 % stats from segmented macrophage mask
%                 curSegPixelCount = sum( imMacrophage(ptCurBinPixelInd) > 0 );
%                 
%                 strField = sprintf('macrophage.segPixelCountOutsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount;', strField));
% 
%                 strField = sprintf('macrophage.segPixelVolOutsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount * volVoxel;', strField));
% 
%                 strField = sprintf('macrophage.segPixelCount.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount + cellColocStats(cellId).macrophage.segPixelCountInsideNucleus;', strField));
%                 
%                 strField = sprintf('macrophage.segPixelVolwithin_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = (curSegPixelCount * volVoxel) + cellColocStats(cellId).macrophage.segPixelVolInsideNucleus;', strField));
% 
%                 % intensity stats within ribbons
%                 strField = sprintf('macrophage.intensityStats.outsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imMacrophageAdjusted, ptCurBinPixelInd);', strField));
% 
%                 strField = sprintf('macrophage.intensityStats.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imMacrophageAdjusted, union(ptCurBinPixelInd, curCellPixelIdxList));', strField));
%                 
%             end
%             
%         % measurments from drug channel 
% 
%             % amount of drug inside cell nucleus
%             cellColocStats(cellId).drug.segPixelCountInsideNucleus = sum(imDrugSeg(curCellPixelIdxList) > 0);
%             cellColocStats(cellId).drug.segPixelVolInsideNucleus = sum(imDrugSeg(curCellPixelIdxList) > 0) * volVoxel;
%             cellColocStats(cellId).drug.segPixelFractionInsideNucleus = mean(imDrugSeg(curCellPixelIdxList) > 0);
% 
%             % drug intensity stats
%             cellColocStats(cellId).drug.intensityStats.insideNucleus = computeIntensityStats(imDrugAdjusted, curCellPixelIdxList);
%             
%             flagInValidNeighborhood = true(numel(imDrug), 1);
%             flagInValidNeighborhood(curCellPixelIdxList) = false;
%             for j = 1:ndims(imMacrophage)
%                 flagInValidNeighborhood(ptPixelLocPhysp(:,j) <= ptCellMinCoord(j) - PARAMETERS.maxDrugNeighDist | ...
%                                         ptPixelLocPhysp(:,j) >= ptCellMaxCoord(j) + PARAMETERS.maxDrugNeighDist) = false;
%             end
%             ptCurPixelLocPhysp = ptPixelLocPhysp(flagInValidNeighborhood, :);
%             ptCurPixelInd = ptPixelInd(flagInValidNeighborhood);
%             
%             [ptIdNearest, dist] = knnsearch(kdtree_ptcell, ptCurPixelLocPhysp);
%             flagValidDist = (dist > 0 & dist < PARAMETERS.maxDrugNeighDist);
%             ptCurPixelInd = ptCurPixelInd(flagValidDist);
%             dist = dist(flagValidDist);
%             
%             binEdges = linspace(0, PARAMETERS.maxDrugNeighDist, PARAMETERS.numDrugNeighDistLevels+1);
%             
%             for bid = 1:numel(binEdges)-1
% 
%                 flagCurBin = (dist >= binEdges(bid) & dist < binEdges(bid+1));
%                 ptCurBinPixelInd = ptCurPixelInd(flagCurBin);
% 
%                 % stats from segmented drug mask
%                 curSegPixelCount = sum( imDrug(ptCurBinPixelInd) > 0 );
%                 
%                 strField = sprintf('drug.segPixelCountOutsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount;', strField));
% 
%                 strField = sprintf('drug.segPixelVolOutsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount * volVoxel;', strField));
% 
%                 strField = sprintf('drug.segPixelCount.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = curSegPixelCount + cellColocStats(cellId).drug.segPixelCountInsideNucleus;', strField));
%                 
%                 strField = sprintf('drug.segPixelVolwithin_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = (curSegPixelCount * volVoxel) + cellColocStats(cellId).drug.segPixelVolInsideNucleus;', strField));
% 
%                 % intensity stats within ribbons
%                 strField = sprintf('drug.intensityStats.outsideNucleus.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imDrugAdjusted, ptCurBinPixelInd);', strField));
% 
%                 strField = sprintf('drug.intensityStats.within_%d_um', binEdges(bid+1));
%                 eval(sprintf('cellColocStats(cellId).%s = computeIntensityStats(imDrugAdjusted, union(ptCurBinPixelInd, curCellPixelIdxList));', strField));
%                 
%             end            
%         
%         % update progress/status
%         percent_done = round(100*cellId/numCells);       
%         
%         if percent_done > last_percent_done
%             fprintf( '%.2d%%  ', percent_done );
%             last_percent_done = percent_done;
%             numPrint = numPrint + 1;
%             if mod( numPrint, 10 ) == 0
%                fprintf( '\n' ); 
%             end
%         end        
%         
%     end
                                            
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
                                                