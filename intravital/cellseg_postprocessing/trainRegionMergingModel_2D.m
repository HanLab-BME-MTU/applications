clc
clear 
close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    defaultAnnotationDataDir = 'C:\deepak\data\regionMergingTraining2D\annotation';
    
    defaultOutputDir = 'C:\deepak\data\regionMergingTraining2D\features';    
    
    % less-critical parameters (leave defaults)
    PARAMETERS.numHistQuantLevels = 64;
    PARAMETERS.bndContStrelRad = 5;
    
%***********************************************************************

% get annotation files
annotationFileList = uipickfiles( 'FilterSpec', fullfile( defaultAnnotationDataDir, '*.rann' ) );

% get output dir
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% generate features for each file
featureData = [];
featureClass = [];
featureStruct = [];

for fid = 1:numel(annotationFileList)

    fprintf( '\n\n*****************************************************\n\n' );           
    
    fprintf( '\n\nComputing Features for file %d/%d ...\n\n', fid, numel(annotationFileList) );
    
    curAnnotationFile = annotationFileList{fid}
    
    % load data from annotation file
    curAnnotationData = load(curAnnotationFile, '-mat', 'annotationData' );    
    curAnnotationData = curAnnotationData.annotationData;
    
    imInput = curAnnotationData.imInput;
    imThresh = curAnnotationData.imThresh;
    imCellSeedPoints = curAnnotationData.imCellSeedPoints;
    imLabelCellSeg = curAnnotationData.imLabelCellSeg;
    numCells = curAnnotationData.numCells
    adjMatrix = curAnnotationData.adjMatrix;
    adjMatrixMerge = curAnnotationData.adjMatrixMerge;
    
    numEdges = 0.5 * numel( find( adjMatrix ) )
    
    % standardize the intensity range of the image 
    imageDynamicRange = ComputeImageDynamicRange( imInput, 99.9 );
    imAdjusted = mat2gray(imInput, imageDynamicRange);    

    % apply median filter to remove any spotty noise
    imAdjusted = medfilt2( imAdjusted, [3,3] );    
    
    % compute basic properties of cell region
    cellStats = [];    
    for cellId = 1:numCells                       

        regProps = ComputeRegionProperties(imLabelCellSeg, cellId);

        % add seed point location
        indCurRegionSeedPoint = regProps.PixelIdxList( imCellSeedPoints(regProps.PixelIdxList) > 0 );
        if isempty( indCurRegionSeedPoint )            
            regProps.ptCellSeedPoint = regProps.Centroid;
        else
            [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
            indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

            regProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
        end

        cellStats = [cellStats; regProps];

    end
    
    % compute region merging features
    [v1, v2, ~] = find( adjMatrix );   
    
    curFeatureData = [];
    curFeatureClass = [];
    
    for k = (find(v1 < v2))'
       
        c1 = v1(k);
        c2 = v2(k);
        
        curFeatureDataStruct = ComputeRegionMergingFeatures(imAdjusted, imLabelCellSeg, cellStats, c1, c2);          
        [curFeatureVec , featureNameList] = ConvertFeatureStructToFeatureVec( curFeatureDataStruct );        
        
        curFeatureData = cat(1, curFeatureData, curFeatureVec');
        
        if adjMatrixMerge(c1,c2) > 0
            curClass = { 'Merge' };
        else
            curClass = { 'Dont-Merge' };
        end
        
        curFeatureClass = cat(1, curFeatureClass, curClass );
        
        curFeatureStruct.annotationFile = curAnnotationFile;
        curFeatureStruct.edgeNodeInd = [c1, c2];
        curFeatureStruct.featureData = curFeatureDataStruct;
        curFeatureStruct.classLabel = adjMatrixMerge(c1,c2);
        curFeatureStruct.featureVec = curFeatureVec;
        curFeatureStruct.featureNameList = featureNameList;
        
        featureStruct = cat(2, featureStruct, curFeatureStruct);
        
    end
    
    fprintf('\nClass Distribution: \n' );
    tabulate( curFeatureClass )
    
    featureData = cat(1, featureData, curFeatureData);
    featureClass = cat(1, featureClass, curFeatureClass);
    
end

fprintf( '\n\n*****************************************************\n\n' );           

fprintf('\nTotal Class Distribution (1 = Merge, -1 = Dont Merge)\n' );
tabulate( featureClass )

numFeatures = size(featureData, 2)
featureNameList

% write csv file
[outFileName, outFilePath] = uiputfile( fullfile(outputRootDir, 'regionMergingFeatures.csv') );
featureLabels = cat(1, 'ClassLabel', featureNameList );
featureMatrix = cat(2, featureClass, featureData );
WriteFeatureMatrixToCSVFile( fullfile(outFilePath, outFileName), featureMatrix, featureLabels );

% write mat file
save( fullfile(outFilePath, strrep( outFileName, '.csv', '.mat' )), 'featureStruct', 'annotationFileList' );

% write arff file for weka
if isempty( cell2mat( strfind( javaclasspath('-dynamic'), 'weka.jar' ) ) )
    javaaddpath( fullfile('C:\Program Files\Weka-3-7', 'weka.jar') );
end
ConvertCSVFileToArffFile( fullfile(outFilePath, outFileName) );

% Train SVM classifier and save the model
% cvp = cvpartition(featureClass, 'kfold', 10);
% 
% svmcrossfun = @(xTrain, yTrain, xTest, rbfSigma, boxConstraint) ...
%                 (svmclassify( ...
%                     svmtrain(xTrain, yTrain, 'Kernel_Function', 'rbf', 'rbf_sigma', rbfSigma, 'boxconstraint', boxConstraint), ...
%                     xTest));
%                 
% svmParamOptFunc = @(z) crossval('mcr', cell2mat(featureData), featureClass, 'partition', cvp, ...
%                                 'Predfun', @(xTrain, yTrain, xTest) svmcrossfun(xTrain, yTrain, xTest, z(1), z(2)) );
% 
%                             
% %[rbfSigmaGrid, boxConstraintGrid] = meshgrid(-15:2:3, -5:2:15);                             
% [rbfSigmaGrid, boxConstraintGrid] = meshgrid(-5:5, -3:2:7);                             
% errorValGrid = zeros(size(rbfSigmaGrid));       
% numParamTrials = numel(rbfSigmaGrid);
% 
% parfor pid = 1:numParamTrials
%     
%     curRbfSigma = rbfSigmaGrid(pid);
%     curBoxConstraint = boxConstraintGrid(pid);    
%     errorValGrid(pid) = svmParamOptFunc( [2^curRbfSigma, 2^curBoxConstraint] );
% 
%     fprintf( '\nGrid Search Trail %d/%d (%.2f%%): Gamma = 2^%d, C = 2^%d, error = %f', ...
%               pid, numParamTrials, 100 * pid / numParamTrials, ...
%               curRbfSigma, curBoxConstraint, errorValGrid(pid));
%           
% end
% 
% [minError, minInd] = min( errorValGrid(:) );
% 
% fprintf( '\nBest Parameters: Gamma = 2^%d, C = 2^%d, error = %f\n', ...
%             rbfSigmaGrid(minInd), boxConstraintGrid(minInd), minError );
%         
% bestRbfSigma = 2^rbfSigmaGrid(minInd)
% bestBoxConstraint = 2^boxConstraintGrid(minInd)
% 
% svmModel = svmtrain( cell2mat(featureData), featureClass, 'Kernel_Function', 'rbf', 'rbf_sigma', bestRbfSigma, 'boxconstraint', bestBoxConstraint);
% save( fullfile(outFilePath, strrep( outFileName, '.csv', '.rmod' )), 'svmModel' );
