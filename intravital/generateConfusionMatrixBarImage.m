
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    validationFileDefaultDir = 'Z:\intravital\data\image_analysis\final_validation';
    
    imageHeightInCells = 85;
    randSeed = 1;
    selFrac = 0.25;

    gapCategory = 50;
    gapLegend = 50;
    widthLegend = 100;
    
    classNameList = { 'G1', 'S', 'G2', 'M' };
    classColors = [235, 64, 61; 242, 234, 73; 117 190 86; 38 169 224];
    
    TrueClassLabel_ColumnName = 'TrueClassLabel';
    PredictedClassLabel_ColumnName = 'PredictedClassLabel';
    PredictionProbability_ColumnName = 'predictionProbability';
    histoneImage_ColumnName = 'histoneImageRelPath';
    fucciImage_ColumnName = 'fucciImageRelPath';
    mipImage_ColumnName = 'mipImageRelPath';
    
    outFileFormat = 'png';
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

validationFileDir = uigetdir( validationFileDefaultDir );

if isempty(validationFileDir)
    return;
end

resultsDir = uigetdir( validationFileDir );

if isempty(resultsDir)
    return;
end

[numData, txtData, rawData] = xlsread( fullfile(validationFileDir, 'cellCycleStateClassificationPerformance.xlsx') );
validationFileHeader = rawData(1,:);
rawData = rawData(2:end,:);

numCells = size(rawData, 1)

trueLabelColumnId = find( strcmpi(validationFileHeader, TrueClassLabel_ColumnName) ) 
predictedLabelColumnId = find( strcmpi(validationFileHeader, PredictedClassLabel_ColumnName) ) 
predictionProbabilityColumnId = find( strcmpi(validationFileHeader, PredictionProbability_ColumnName) ) 
histoneImageColumnId = find( strcmpi(validationFileHeader, histoneImage_ColumnName) ) 
fucciImageColumnId = find( strcmpi(validationFileHeader, fucciImage_ColumnName) ) 
mipImageColumnId = find( strcmpi(validationFileHeader, mipImage_ColumnName) ) 

% compute confusion matrix
cMatrix = zeros( numel(classNameList), numel(classNameList) );
cMatrixNmzd = zeros( numel(classNameList), numel(classNameList) );
cMatrixCellList = cell(numel(classNameList), numel(classNameList));
incCellList = [];

for cidTrue = 1:numel(classNameList)
    
    curTrueClass = classNameList{cidTrue};
    curTrueCellList = find( strcmp(rawData(:,trueLabelColumnId), curTrueClass) );
    
    for cidPred = 1:numel(classNameList)
        
        curPredClass = classNameList{cidPred};        
        curPredCellList = curTrueCellList( strcmp(rawData(curTrueCellList,predictedLabelColumnId), curPredClass) );
        
        cMatrix(cidTrue, cidPred) = numel(curPredCellList);
        cMatrixNmzd(cidTrue, cidPred) = numel(curPredCellList) / numel(curTrueCellList);

    end
    
end

cMatrix
cMatrixNmzd

% find number of rows and columns for each category in the confusion matrix
if selFrac < 1
    rgen = RandStream('mt19937ar','Seed',randSeed);
    rgen.reset();
    selCellList = randsample(rgen, numCells, round(selFrac * numCells));
else
    selCellList = 1:numCells;
end

barImageSizeInCells = zeros(numel(classNameList), numel(classNameList), 2);

for cidTrue = 1:numel(classNameList)
    
    curTrueClass = classNameList{cidTrue};
    
    curTrueCellList = selCellList( strcmp(rawData(selCellList,trueLabelColumnId), curTrueClass) );
    numCellsInCurClass = numel(curTrueCellList)
    
    maxColumns = 0;
    
    for cidPred = 1:numel(classNameList)
        
        curPredClass = classNameList{cidPred};        
        curPredCellList = curTrueCellList( strcmp(rawData(curTrueCellList,predictedLabelColumnId), curPredClass) );
        
        numCellsInColumn = ceil(imageHeightInCells * cMatrixNmzd(cidTrue, cidPred));     
        numColumns = ceil( numel(curPredCellList)/numCellsInColumn );
        
        if numColumns > maxColumns 
            maxColumns = numColumns
        end
        
        barImageSizeInCells(cidTrue, cidPred, 1) = numCellsInColumn;
        
    end

    for cidPred = 1:numel(classNameList)        
        barImageSizeInCells(cidTrue, cidPred, 2) = maxColumns;    
    end

    if sum( barImageSizeInCells(cidTrue, :, 1) ) > imageHeightInCells        
        [maxval, maxCatInd] = max(barImageSizeInCells(cidTrue, :, 1));
        diffRows = sum( barImageSizeInCells(cidTrue, :, 1) ) - imageHeightInCells;
        barImageSizeInCells(cidTrue, maxCatInd, 1) = barImageSizeInCells(cidTrue, maxCatInd, 1) - diffRows;
    end
    
end

% generate images for each cell class
thumbImageColIdList = [histoneImageColumnId, fucciImageColumnId, mipImageColumnId];

thumbsize = size( imread(fullfile(validationFileDir, rawData{1, histoneImageColumnId})) );

numTotalCellsDisplayed = 0;

for cidTrue = 1:numel(classNameList)

    curTrueClass = classNameList{cidTrue};
    
    fprintf( '\n>> %d/%d: Generating image for %s\n', cidTrue, numel(classNameList), curTrueClass);
    
    curTrueCellList = find( strcmp(rawData(:,trueLabelColumnId), curTrueClass) );
    numCellsInCurClass = numel(curTrueCellList)
    
    imCurClassBar = cell(1, numel(thumbImageColIdList));
    imCurLegendBar = [];
    
    for cidPred = 1:numel(classNameList)
    
        curPredClass = classNameList{cidPred};        
        fprintf( '\n\nPredicted as %s: \n', curPredClass );
        curPredCellList = curTrueCellList( strcmp(rawData(curTrueCellList,predictedLabelColumnId), curPredClass) );

        numCellsToSelect = prod( barImageSizeInCells(cidTrue, cidPred, :) );
        
        if numel(curPredCellList) > numCellsToSelect             
            curPredCellList = curPredCellList( randsample(rgen, numel(curPredCellList), numCellsToSelect) );
        end
        
        numCellsSelected = numel(curPredCellList)        
        curBarImageSizeInCells = (squeeze(barImageSizeInCells(cidTrue, cidPred, :)))'
        
        numTotalCellsDisplayed = numTotalCellsDisplayed + numCellsSelected;
        
        for thid = 1:numel(thumbImageColIdList)
            
            imCurThumbnailBar = zeros( [curBarImageSizeInCells .* thumbsize(1,2), 3], 'uint8');  
            r = 0;
            c = 0;
            
            fprintf( '\nProgress: \n' );
            last_percent_done = 0;
            numPrint = 0;
            
            for i = 1:numel(curPredCellList) 
                
                cid = curPredCellList(i);

                imCurThumbnail = imread( fullfile(validationFileDir, rawData{cid, thumbImageColIdList(thid)}) );

                if ndims(imCurThumbnail) < 3
                    imCurThumbnail = repmat(imCurThumbnail, [1, 1, 3]);
                end

                imagePos = {r * thumbsize(1) + (1:thumbsize(1)), c * thumbsize(1) + (1:thumbsize(2)), ':'};
                imCurThumbnailBar( imagePos{:} ) = imCurThumbnail;
                
                r = r+1;
                if r >= curBarImageSizeInCells(1)
                    r = 0;
                    c = c + 1;
                end

                percent_done = round(100*i/numel(curPredCellList));       

                if percent_done > last_percent_done
                    fprintf( '%.2d%%  ', percent_done );
                    last_percent_done = percent_done;
                    numPrint = numPrint + 1;
                    if mod( numPrint, 10 ) == 0
                       fprintf( '\n' ); 
                    end
                end        
                
            end
            
            if cidPred > 1
                imCurThumbnailBar = padarray(imCurThumbnailBar, gapCategory * [1 0 0], 255, 'pre');
            end

            imCurClassBar{thid} = cat(1, imCurClassBar{thid}, imCurThumbnailBar);
            
        end

        imCurLegend = repmat( reshape(classColors(cidPred, :), 1, 1, []), curBarImageSizeInCells(1) * thumbsize(1), widthLegend);
        if cidPred > 1
            imCurLegend = padarray(imCurLegend, gapCategory * [1 0 0], 255, 'pre');
        end
        imCurLegendBar = cat(1, imCurLegendBar, imCurLegend);
        
    end

    % save individual thumbnail bar images
    for thid = 1:numel(thumbImageColIdList)
        curOutFileName = fullfile(resultsDir, sprintf('%s_%d.%s', curTrueClass, thid, outFileFormat));
        imwrite( imCurClassBar{thid}, curOutFileName, outFileFormat );
    end

    % save legend image
    curOutFileName = fullfile(resultsDir, [curTrueClass, '_legend.' outFileFormat]);
    imwrite( uint8(imCurLegendBar), curOutFileName, outFileFormat );
    
    % save bar image
    imCurClassBarCombined = cat(2, imCurClassBar{:});             
    curOutFileName = fullfile(resultsDir, [curTrueClass, '.' outFileFormat]);
    imwrite( imCurClassBarCombined, curOutFileName, outFileFormat );

    % save bar image with legend
    imCurClassBarCombined_with_legend = cat(2, imCurClassBarCombined, padarray(imCurLegendBar, gapLegend * [0, 1, 0], 255, 'pre'));
    curOutFileName = fullfile(resultsDir, [curTrueClass, '_with_legend.' outFileFormat]);
    imwrite( imCurClassBarCombined_with_legend, curOutFileName, outFileFormat );
    
end

fprintf( '\nA total of %d/%d (%.2f%%) cells are displayed ... \n', ...
         numTotalCellsDisplayed, numCells, ...
         100.0 * numTotalCellsDisplayed/numCells);
     
% % compute confusion matrix
% cMatrix = zeros( numel(classNameList), numel(classNameList) );
% cMatrixNmzd = zeros( numel(classNameList), numel(classNameList) );
% cMatrixCellList = cell(numel(classNameList), numel(classNameList));
% incCellList = [];
% 
% for cidTrue = 1:numel(classNameList)
%     
%     curTrueClass = classNameList{cidTrue};
%     curTrueCellList = find( strcmp(rawData(:,trueLabelColumnId), curTrueClass) );
%     
%     for cidPred = 1:numel(classNameList)
%         
%         curPredClass = classNameList{cidPred};        
%         curPredCellList = curTrueCellList( strcmp(rawData(curTrueCellList,predictedLabelColumnId), curPredClass) );
%         
%         cMatrix(cidTrue, cidPred) = numel(curPredCellList);
%         cMatrixNmzd(cidTrue, cidPred) = numel(curPredCellList) / numel(curTrueCellList);
% 
%         if numel(curPredCellList) < minCellCount
%             incCellList = [ incCellList; curPredCellList ];
%         end
%     end
%     
% end
% 
% cMatrix
% cMatrixNmzd
% 
% % choose random subset
% if selFrac < 1
%     rgen = RandStream('mt19937ar','Seed',randSeed);
%     rgen.reset();
%     selCellList = randsample(rgen, numCells, round(selFrac * numCells));
%     selCellList = union(selCellList, incCellList);
% else
%     selCellList = 1:numCells;
% end
% 
% % generate images for each cell class
% thumbImageColIdList = [histoneImageColumnId, fucciImageColumnId, mipImageColumnId];
% 
% thumbsize = size( imread(fullfile(validationFileDir, rawData{1, histoneImageColumnId})) );
% 
% for cidTrue = 1:numel(classNameList)
% 
%     curTrueClass = classNameList{cidTrue};
%     
%     fprintf( '\n>> %d/%d: Generating image for %s\n', cidTrue, numel(classNameList), curTrueClass);
%     
%     curTrueCellList = selCellList( strcmp(rawData(selCellList,trueLabelColumnId), curTrueClass) );
%     numCellsInCurClass = numel(curTrueCellList)
%     
%     imCurClassBar = cell(numel(classNameList), numel(thumbImageColIdList));
%     maxColumns = 0;
%     
%     for cidPred = 1:numel(classNameList)
%     
%         curPredClass = classNameList{cidPred};        
%         fprintf( '\nPredicted as %s: \n', curPredClass );
%         curPredCellList = curTrueCellList( strcmp(rawData(curTrueCellList,predictedLabelColumnId), curPredClass) );
% 
%         %numTotalCellsInColumn = ceil(imageHeightInCells * cMatrix(cidTrue, cidPred))
%         numTotalCellsInColumn = round(imageHeightInCells * numel(curPredCellList) / numel(curTrueCellList))
%         if numTotalCellsInColumn < 1
%             numTotalCellsInColumn = 1;
%         end
%         
%         numColumns = ceil( numel(curPredCellList)/numTotalCellsInColumn )         
%         
%         for thid = 1:numel(thumbImageColIdList)
%             
%             imCurThumbnailBar = zeros(numTotalCellsInColumn * thumbsize(1), numColumns * thumbsize(2), 3, 'uint8');  
%             r = 0;
%             c = 0;
%             
%             fprintf( '\nProgress: \n' );
%             last_percent_done = 0;
%             numPrint = 0;
%             
%             for i = 1:numel(curPredCellList) 
%                 
%                 cid = curPredCellList(i);
% 
%                 imCurThumbnail = imread( fullfile(validationFileDir, rawData{cid, thumbImageColIdList(thid)}) );
% 
%                 if ndims(imCurThumbnail) < 3
%                     imCurThumbnail = repmat(imCurThumbnail, [1, 1, 3]);
%                 end
% 
%                 imagePos = {r * thumbsize(1) + (1:thumbsize(1)), c * thumbsize(1) + (1:thumbsize(2)), ':'};
%                 imCurThumbnailBar( imagePos{:} ) = imCurThumbnail;
%                 
%                 r = r+1;
%                 if r >= numTotalCellsInColumn
%                     r = 0;
%                     c = c + 1;
%                 end
% 
%                 percent_done = round(100*i/numel(curPredCellList));       
% 
%                 if percent_done > last_percent_done
%                     fprintf( '%.2d%%  ', percent_done );
%                     last_percent_done = percent_done;
%                     numPrint = numPrint + 1;
%                     if mod( numPrint, 10 ) == 0
%                        fprintf( '\n' ); 
%                     end
%                 end        
%                 
%             end
%             
%             imCurClassBar{cidPred, thid} = imCurThumbnailBar;
%             
%         end
%         
%         if size(imCurClassBar{cidPred, 1}, 2) > maxColumns
%             maxColumns = size(imCurClassBar{cidPred, 1}, 2);
%         end
%         
%     end
% 
%     imCurClassBarCombined = cell(1, numel(thumbImageColIdList));
%     for thid = 1:numel(thumbImageColIdList)
%        
%         for cidPred = 1:numel(classNameList)
%             imCurClassBar{cidPred, thid}(:,end+1:maxColumns,:) = 0;
%             imCurClassBar{cidPred, thid} = padarray(imCurClassBar{cidPred, thid}, gapCategory * [1 0 0], 255);
%         end
%         
%         imCurClassBarCombined{thid} = cat(1, imCurClassBar{:,thid});
%         
%         curOutFileName = fullfile(resultsDir, sprintf('%s_%d.%s', curTrueClass, thid, outFileFormat));
%         imwrite( imCurClassBarCombined{thid}, curOutFileName, outFileFormat );
%         
%     end
%     
%     imCurClassBarCombined = cat(2, imCurClassBarCombined{:});         
%     
%     curOutFileName = fullfile(resultsDir, [curTrueClass, '.' outFileFormat]);
%     imwrite( imCurClassBarCombined, curOutFileName, outFileFormat );
%     
% end

