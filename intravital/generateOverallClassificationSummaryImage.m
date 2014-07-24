
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    validationFileDefaultDir = 'Z:\intravital\data\image_analysis\final_validation';
    
    numProbBins = 4;
    numImagesPerBin = 6;

    minProb = -1;
    selModeList = {'random', 'decprob'};
    selMode = 1;
    randSeed = 6; % good values 0.8 --- 5, 6; -1 -- 6
    
    probBoxWidthFrac = 0.6;
    probBoxBorderWidth = 6.0;
    probBoxGapSize = 4.0;
    probPadSize = 4;
    imageGapSize = 0;
    
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

classNameList = unique(rawData(:,predictedLabelColumnId));
predictionProbabilities = cell2mat(rawData(:, predictionProbabilityColumnId));
predictionProbabilityDist = cell2mat(rawData(:, predictionProbabilityColumnId+1:predictionProbabilityColumnId+numel(classNameList)));

% generate images for each cell cycle state
rgen = RandStream('mt19937ar','Seed',randSeed);
rgen.reset();

for i = 1:numel(classNameList)
   
    curClass = classNameList{i};
    
    fprintf( '\n%d/%d: Generating image for %s\n', i, numel(classNameList), curClass);
    
    % find all cells in the current category
    curClassCellList = find( strcmp(rawData(:,predictedLabelColumnId), curClass) );

    % sort cells in decreasing order of probability
    curClassCellList = sortrows([curClassCellList,  -1 * predictionProbabilities(curClassCellList)], 2);
    curClassCellList = curClassCellList(:, 1);

    % bin cells by probability and generate an image for each bin
    curPredictionProbabilities = predictionProbabilities(curClassCellList);

    binWidth = 1 / numProbBins;
    
    if minProb > 0
        probRange = [minProb, max(curPredictionProbabilities)];
    else
        probRange = [min(curPredictionProbabilities), max(curPredictionProbabilities)];
    end
    
    curPredictionProbabilityBins = floor( mat2gray(curPredictionProbabilities, probRange) / binWidth ) + 1;
    curPredictionProbabilityBins(curPredictionProbabilityBins > numProbBins) = numProbBins;
    
    imCurClassSummary = [];
    
    for bid = 1:numProbBins             
    
        curBinCellList = curClassCellList( curPredictionProbabilityBins == bid );
        
        % select a subset of cells of the specified size
        switch selModeList{selMode}

            case 'random'

                selCellList = curBinCellList(randperm(rgen, numel(curBinCellList), min(numImagesPerBin, numel(curBinCellList))));
                selCellList = sortrows([selCellList,  -1 * predictionProbabilities(selCellList)], 2);
                selCellList = selCellList(:, 1);
                
            case 'decprob'

                selCellList = curBinCellList(1:min(numImagesPerBin, numel(curBinCellList)));
                
        end
        
        % generate image from the selected cells
        imCurBinSummary = [];

        for j = 1:numel(selCellList)

            curCellId = selCellList(j); 

            imHistone = imread( fullfile(validationFileDir, rawData{curCellId, histoneImageColumnId}) );
            imFucci = imread( fullfile(validationFileDir, rawData{curCellId, fucciImageColumnId}) );
            imMIP = imread( fullfile(validationFileDir, rawData{curCellId, mipImageColumnId}) );
            imMIP = repmat(imMIP, [1, 1, 3]);

            imCurCellSummary = cat(1, imHistone, imFucci, imMIP);

            % shown probability distribution
            imProb = [];
            probBoxSize = floor(size(imHistone,2) * probBoxWidthFrac * 0.25);
            curClassProbabilityDist = predictionProbabilityDist(curCellId, :);
            for lid = 1:numel(classNameList)
               imCurProb = repmat( (1-curClassProbabilityDist(lid)) * 255 * ones( probBoxSize * [1,1] ), [1, 1, 3]);
               imCurProb = padarray(imCurProb, probBoxBorderWidth * [1, 1, 0], 0.0);
               imCurProb = padarray(imCurProb, probBoxGapSize * [1, 1, 0], 255.0);
               imProb = cat(2, imProb, imCurProb);
            end
            imProb = padarray(imProb, probPadSize * [1, 1, 0], 255.0);
            imProb(:,end:size(imHistone,2),:) = 255;
            imProb = padarray(imProb, 4 * [2, 1, 0], 0.0);
            
            probimsize = size(imProb);
            probscalefactor = size(imHistone,2)/size(imProb,2);
            probimoutsize = round(probscalefactor * probimsize(1:2));
            imProb = imresize(imProb, probimoutsize);
            
            imCurCellSummary = cat(1, imProb, imCurCellSummary);
            
            % add gap between images
            imCurCellSummary = padarray(imCurCellSummary, [imageGapSize, imageGapSize, 0], 255);
            
            % add to strip
            imCurBinSummary = cat(2, imCurBinSummary, imCurCellSummary);

        end        

        imCurClassSummary = cat(2, imCurBinSummary, imCurClassSummary);
        
    end

    curOutFileName = fullfile(resultsDir, [curClass, '.' outFileFormat]);
    imwrite( imCurClassSummary, curOutFileName, outFileFormat );
    
end

