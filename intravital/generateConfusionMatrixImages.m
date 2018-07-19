
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    validationFileDefaultDir = 'Z:\intravital\data\image_analysis\final_validation';
    numImagesPerGroup = 10;
    gapSize = 0;
    
    selModeList = {'random', 'confidence'};
    selMode = 2;
    
    TrueClassLabel_ColumnName = 'TrueClassLabel';
    PredictedClassLabel_ColumnName = 'PredictedClassLabel';
    PredictionProbability_ColumnName = 'predictionProbability';
    histoneImage_ColumnName = 'histoneImageRelPath';
    fucciImage_ColumnName = 'fucciImageRelPath';
    mipImage_ColumnName = 'mipImageRelPath';
    
    outFileFormat = 'tif';
    
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

% compute the confusion matrix category of each cell
cellConfusionMatrixCategory = cell(numCells, 1);
for cid = 1:numCells
    curTrueLabel = rawData{cid, trueLabelColumnId};
    curPredictedLabel = rawData{cid, predictedLabelColumnId};
    cellConfusionMatrixCategory{cid} = [curTrueLabel '_' curPredictedLabel];
end

confusionMatrixCategories = unique(cellConfusionMatrixCategory)

% generate images for each category of the confusion matrix
for i = 1:numel(confusionMatrixCategories)
   
    curCategory = confusionMatrixCategories{i};
    
    fprintf( '\n%d/%d: Generating image for %s\n', i, numel(confusionMatrixCategories), curCategory);
    
    % find all cells in the current category
    curCellList = find( strcmp(cellConfusionMatrixCategory, curCategory) );

    % select a specified number of cells randomly
    switch selModeList{selMode}

        case 'random'
            
            selCellList = curCellList(randperm(numel(curCellList), min(numImagesPerGroup, numel(curCellList))));
            
        case 'confidence'
            
            selCellList = sortrows([curCellList,  -1 * cell2mat(rawData(curCellList, predictionProbabilityColumnId))], 2);
            selCellList = selCellList(1:min(numImagesPerGroup, numel(curCellList)), 1 );

    end
    
    % generate image from the selected cells
    imCurCategory = [];
    
    for j = 1:numel(selCellList)
        
        curCellId = selCellList(j); 

        imHistone = imread( fullfile(validationFileDir, rawData{curCellId, histoneImageColumnId}) );
        imFucci = imread( fullfile(validationFileDir, rawData{curCellId, fucciImageColumnId}) );
        imMIP = imread( fullfile(validationFileDir, rawData{curCellId, mipImageColumnId}) );
        imMIP = repmat(imMIP, [1, 1, 3]);
        
        imCellSummary = padarray(cat(1, imHistone, imFucci, imMIP), [gapSize, gapSize, 0], 255);
        
        imCurCategory = cat(2, imCurCategory, imCellSummary);
        
    end
    
    curOutFileName = fullfile(resultsDir, [curCategory, '.' outFileFormat]);
    imwrite( imCurCategory, curOutFileName, outFileFormat );
    
end

