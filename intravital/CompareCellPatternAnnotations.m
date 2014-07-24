function [performanceStats] = CompareCellPatternAnnotations( annotationFile1, annotationFile2 )

    if ~exist( 'annotationFile1', 'var' )       
        [fileName,pathName] = uigetfile( fullfile(pwd, 'CellPatternAnnotation.mat' ), 'Select annotation file - 1' );   
        annotationFile1 = fullfile(pathName, fileName);
    end
    
    if ~exist( 'annotationFile2', 'var' )       
        [fileName,pathName] = uigetfile( fullfile(pwd, 'CellPatternAnnotation.mat' ), 'Select annotation file - 1' );   
        annotationFile2 = fullfile(pathName, fileName);
    end
    
    annotationData1 = load( annotationFile1 );
    annotationData2 = load( annotationFile2 );

    assert( numel(annotationData1.cellStats) == numel(annotationData2.cellStats) );
    
    numCells = numel( annotationData1.cellStats );
    
    cellClass1 = cell( numCells, 1 );
    cellClass2 = cell( numCells, 1 );
    flagIsValidCell = true(numCells, 1);
    
    for i = 1:numCells
        
        cellClass1{i} = DetermineCellClass( annotationData1.cellStats(i).cellPatternType );
        cellClass2{i} = DetermineCellClass( annotationData2.cellStats(i).cellPatternType );
        
        if isempty(cellClass1{i}) || isempty(cellClass2{i})
            flagIsValidCell(i) = false;
        end
        
    end
    
    labelList = { 'G1', 'S', 'G2', 'M' };
    performanceStats = ComputeClassificationPerformance_MultiClass( cellClass1(flagIsValidCell), cellClass2(flagIsValidCell), labelList );
    
    performanceStats.cMatrix
    
end

function [cellClass] = DetermineCellClass( cellPatternType )

   switch(cellPatternType)

        case { 'Mono_Pre_G1', 'Multi_Pre_G1', 'Mono_G1', 'Multi_G1' } 

            cellClass = 'G1';

        case { 'Mono_S', 'Multi_S' } 

            cellClass = 'S';

        case { 'Mono_G2', 'Multi_G2' }     

            cellClass = 'G2';

        case { 'Mono_Prophase', 'Multi_Prophase' }     

            cellClass = 'M';                    

        otherwise

            cellClass = [];
            
   end             


end



