function [cellColocStats] = ComputeDNADamageColocalizationMeasures(imDrug, im53BP1, imMacrophage, spacing, ...
                                                                   imLabelCellSeg, cellStats, ...
                                                                   imDrugSeg, imMacrophageSeg, varargin)
                                            
    p = inputParser;
    p.addParamValue('maxDrugNeighDist', 50, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('numDrugNeighDistLevels', 5, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse(varargin{:});
          
    PARAMETERS = p.Results;
    
    numCells = numel(cellStats);
    volVoxel = prod(spacing);
    
    imMacrophageDistMap = bwdistsc(imMacrophageSeg, spacing);
    
    drugPixelInd = find(imDrugSeg > 0);
    ptDrugPixelLocPhysp = ind2submat(size(imDrugSeg), drugPixelInd) * diag(spacing);
    
    imLabelCellBoundarySeg = bwperim(imLabelCellSeg > 0) .* imLabelCellSeg;
    
    fprintf( '\nProgress: \n' );
    
    last_percent_done = 0;
    numPrint = 0;
    
    for cellId = 1:numCells
        
        curCellPixelIdxList = cellStats(cellId).PixelIdxList;
        
        cellColocStats(cellId).cellId = cellId;
        
        % distance to nearest macrophage
        cellColocStats(cellId).distanceToMacrophage = min( imMacrophageDistMap(curCellPixelIdxList) );
        
        % amount of drug inside cell nucleus
        cellColocStats(cellId).drugVolInsideNucleus = sum(imDrugSeg(curCellPixelIdxList) > 0) * volVoxel;
        cellColocStats(cellId).drugFractionInsideNucleus = mean(imDrugSeg(curCellPixelIdxList) > 0);
        
        % amount of drug around cell nucleus
        curCellBoundaryPixelIdxList = find(imLabelCellBoundarySeg == cellId);
        ptCurCellBoundaryPixelLocPhysp = ind2submat(size(imLabelCellSeg), curCellBoundaryPixelIdxList) * diag(spacing);
        kdtree_ptcell = KDTreeSearcher(ptCurCellBoundaryPixelLocPhysp);          

        flagInValidNeighborhood = true(numel(drugPixelInd), 1);
        flagInValidNeighborhood(imLabelCellSeg(drugPixelInd) == cellId) = false;
        ptCellMinCoord = min(ptCurCellBoundaryPixelLocPhysp);
        ptCellMaxCoord = max(ptCurCellBoundaryPixelLocPhysp);
        for j = 1:ndims(imDrug)
            flagInValidNeighborhood(ptDrugPixelLocPhysp(:,j) <= ptCellMinCoord(j) - PARAMETERS.maxDrugNeighDist | ...
                                    ptDrugPixelLocPhysp(:,j) >= ptCellMaxCoord(j) + PARAMETERS.maxDrugNeighDist) = false;
        end
        
        [ptIdNearest, dist] = knnsearch(kdtree_ptcell, ptDrugPixelLocPhysp(flagInValidNeighborhood, :));
        dist(dist == 0 | dist >= PARAMETERS.maxDrugNeighDist) = [];
        binEdges = linspace(0, PARAMETERS.maxDrugNeighDist, PARAMETERS.numDrugNeighDistLevels+1);
        binDrugCounts = cumsum(histc(dist, binEdges));
        
        for bid = 1:PARAMETERS.numDrugNeighDistLevels

            strField = sprintf('drugVolOutsideNucleusWithin_%d_um', binEdges(bid+1));
            cellColocStats(cellId).(strField) = binDrugCounts(bid) * volVoxel;
            
            strField = sprintf('drugVolWithin_%d_um', binEdges(bid+1));
            cellColocStats(cellId).(strField) = (binDrugCounts(bid) * volVoxel) + cellColocStats(cellId).drugVolInsideNucleus;
            
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
                                            
                                                