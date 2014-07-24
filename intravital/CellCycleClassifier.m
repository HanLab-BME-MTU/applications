classdef CellCycleClassifier < handle    
    
    properties (SetAccess = private)
        flagModelReady
        G1_S_G2M_Model;
        G2_M_Model;
    end
    
    methods       
        
        function obj = CellCycleClassifier()          

            obj.flagModelReady = false;
            obj.G1_S_G2M_Model = [];
            obj.G2_M_Model = [];
            
            AddWekaClassesToPath();
            
        end
            
        function [flagLoadStatus] = LoadModelFromDisk(obj, modelDir)
            
            % weka imports
            import weka.core.*;
            import weka.core.converters.*;
            import weka.classifiers.*;
            
            assert( exist(fullfile(modelDir, 'G1_S_G2M.model'), 'file') > 0 );
            assert( exist(fullfile(modelDir, 'G2_M.model'), 'file') > 0 );
            
            % load G1_S_G2M model
            modelFileContents = SerializationHelper.readAll( fullfile(modelDir, 'G1_S_G2M.model') );
            obj.G1_S_G2M_Model.classifier =  modelFileContents(1);
            obj.G1_S_G2M_Model.header =  modelFileContents(2);

            % load G2_M model
            modelFileContents = SerializationHelper.readAll( fullfile(modelDir, 'G2_M.model') );
            obj.G2_M_Model.classifier = modelFileContents(1);
            obj.G2_M_Model.header = modelFileContents(2);
            
            obj.flagModelReady = true;
            
        end
        
        function BuildModel(obj, flagIgnoreCellsOnBorder, defaultDataRootDir)            
            
            % weka imports
            import weka.core.*;
            import weka.core.converters.*;

            import weka.classifiers.meta.*;
            import weka.filters.supervised.instance.*;
            import weka.classifiers.functions.*;              
            import weka.classifiers.trees.*;       
            import weka.classifiers.*;

            import java.util.*;
            
            % Get arguments
            if ~exist( 'flagIgnoreCellsOnBorder', 'var' )
                flagIgnoreCellsOnBorder = true;                
            end
            
            if ~exist( 'defaultDataRootDir', 'var' )
                defaultDataRootDir = pwd;
            end
            
            % Ask the user to provide the list of directories to be processed
            rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);

            if isempty(rootDirList)
                return;
            end
            
            % Ask the user to select the output directory
            outputRootDir = uigetdir( defaultDataRootDir, 'Select Model Output Directory' ); 
            
            if isempty(outputRootDir)
                error( 'Model output directory must be specified' );
                return;
            end
            
            % generate features            
            classificationTask = { 'G1_S_G2M', 'G2_M' };

            globalFeatureData = cell(size(classificationTask));
            hWaitBar = waitbar(0, 'Computing Features');
            
            for rid = 1:numel(rootDirList)

                fprintf( '\n\n*****************************************************\n\n' );           

                % get a directory to process from the list
                dataRootDir = rootDirList{rid};     

                fprintf( '\nGenerating features for root directory - %d/%d\n', rid, numel(rootDirList) );

                % set output directory accordingly
                [pathstr, name, ext] = fileparts( dataRootDir );
                outputDir = fullfile(outputRootDir, name);

                if ~isdir(outputDir)
                    mkdir(outputDir);
                end

                dataRootDir
                outputDir

                % log
                diary off;
                diary_file = fullfile( outputDir , sprintf( '%s.log' , mfilename ) );
                if exist( diary_file, 'file' )
                    fclose all;
                    delete( diary_file );
                end
                diary( diary_file );
                diary on;

                % get a list of all cell pattern annotation files
                annotationFileList = rdir( fullfile(dataRootDir, '**/CellPatternAnnotation.mat') );

                % generate features for each annotation file
                for fid = 1:numel(annotationFileList)

                   fprintf( '\n\nComputing Features for root directory %d/%d and file %d/%d ...\n\n', ...
                            rid, numel(rootDirList), fid, numel(annotationFileList) );       

                   tic

                   curAnnotationFilePath = annotationFileList(fid).name   
                   curAnnotationRelativeFilePath = annotationFileList(fid).name((numel(dataRootDir)+2):end);

                   [pathstr, name, ext] = fileparts(curAnnotationRelativeFilePath); 
                   curOutputDir = fullfile( outputDir, pathstr );
                   if ~isdir(curOutputDir)
                       mkdir(curOutputDir);
                   end

                   % load annotation data
                   waitbar(0, hWaitBar, ...
                              sprintf('Loading file %d/%d in root directory %d/%d', ...
                                       fid, numel(annotationFileList), rid, numel(rootDirList) ));
                   
                   try
                        annotationData = load( curAnnotationFilePath ); 
                   catch
                       fprintf('\nERROR: could not load annotation file\n' );
                       continue;
                   end

                   % standardize the image data before computing any features
                   imageDataAdjusted = StandardizeIntravitalImageData(annotationData.imageData);

                   % denoise it a bit using a median filter
                   for i = 1:numel(annotationData.imageData)
                       imageDataAdjusted{i} = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted{i} );    
                   end

                   % make a note of all the cells touching the border
                   borderCellIds = setdiff( unique(annotationData.imLabelCellSeg .* ~annotationData.imRegValidMask), 0 );

                   % compute features for each cell
                   cellFeatureData = cell(size(classificationTask));

                   flagIsCellWellSegmented = zeros(numel(annotationData.cellStats),1);
                   flagIsCellOnBorder = zeros(numel(annotationData.cellStats),1);
                   
                   waitbar(0, hWaitBar, ...
                              sprintf('Computing Features for file %d/%d in root directory %d/%d', ...
                                       fid, numel(annotationFileList), rid, numel(rootDirList) ));
                   
                   for i = 1:numel(annotationData.cellStats)

                        curCellStats = annotationData.cellStats(i);

                        fprintf( '\n>> Computing Features for cell %d/%d of type - %s ...\n', i, numel(annotationData.cellStats), curCellStats.cellPatternType );       

                        % ignore the badly segmented cells
                        if ismember( curCellStats.cellPatternType, { 'Bad_Detection', 'Under_Segmented', 'Over_Segmented'} );
                            flagIsCellWellSegmented(i) = false;
                            fprintf( '\nCell was erroneously segmented/detected\n' );
                            continue;            
                        end

                        flagIsCellWellSegmented(i) = true;       

                        % check if cell touches the image border
                        flagIsCellOnBorder(i) = ismember( i, borderCellIds );       

                        if flagIsCellOnBorder(i)
                            fprintf( '\nBorder Cell - Part of the cell is outside the valid foreground region\n' );
                            if flagIgnoreCellsOnBorder
                                continue;
                            end
                        end

                        % note down some cell indentity information 
                        curCellStats.cellId = i;        
                        curCellStats.flagIsOnBorder = flagIsCellOnBorder(i);        
                        curCellStats.annotationFile = curAnnotationFilePath;        
                        imCurCellSegMask = (annotationData.imLabelCellSeg == i);


                        % Compute features for each classification task
                        for tid = 1:numel(classificationTask)
                            featurefunc = str2func( ['ComputeCellPatternFeatures_', ...
                                                      classificationTask{tid}] );
                                                  
                            [classNameList, featureStruct, ...
                            featureVec, featureLabels, ...
                            flagLearnFeature ] = featurefunc( imageDataAdjusted, ...
                                                              imCurCellSegMask, ...
                                                              curCellStats );
                                                      

                            if isempty(featureStruct)
                                continue;
                            end

                            if isempty(cellFeatureData{tid})
                                
                                cellFeatureData{tid}.featureStruct = featureStruct;
                                cellFeatureData{tid}.featureMatrix = featureVec';
                                cellFeatureData{tid}.featureLabels = featureLabels';
                                cellFeatureData{tid}.flagUseAsFeature = flagLearnFeature;
                                cellFeatureData{tid}.classNameList = classNameList;                                
                                cellFeatureData{tid}.classificationTask = classificationTask{tid};
                                
                            else
                                
                                cellFeatureData{tid}.featureStruct = cat(1, cellFeatureData{tid}.featureStruct, featureStruct);
                                cellFeatureData{tid}.featureMatrix = cat(1, cellFeatureData{tid}.featureMatrix, featureVec');
                                
                            end
                            
                            
                            if isempty(globalFeatureData{tid})
                                
                                globalFeatureData{tid} = cellFeatureData{tid};
                                
                                [wekaInstance, ...
                                 wekaAttributeList, ...
                                 classIndex] = CellCycleClassifier.ConstructWekaInstance( featureVec(flagLearnFeature > 0), ...
                                                                                          featureLabels(flagLearnFeature > 0), ...
                                                                                          classNameList);
                                
                                globalFeatureData{tid}.weka.attributeList = wekaAttributeList;                                
                                globalFeatureData{tid}.weka.trainingDataset = weka.core.Instances( [classificationTask{tid} '_Training'], ...
                                                                                                   globalFeatureData{tid}.weka.attributeList, ...
                                                                                                   numel(annotationData.cellStats) );
                                
                                globalFeatureData{tid}.weka.trainingDataset.add( wekaInstance );
                                globalFeatureData{tid}.weka.trainingDataset.setClassIndex(classIndex - 1);
                                
                            else
                                globalFeatureData{tid}.featureStruct = cat(1, globalFeatureData{tid}.featureStruct, cellFeatureData{tid}.featureStruct);
                                globalFeatureData{tid}.featureMatrix = cat(1, globalFeatureData{tid}.featureMatrix, cellFeatureData{tid}.featureMatrix);
                                
                                [wekaInstance, ~, ~] = CellCycleClassifier.ConstructWekaInstance( featureVec(flagLearnFeature > 0), ...
                                                                                                  featureLabels(flagLearnFeature > 0), ...
                                                                                                  classNameList);
                                globalFeatureData{tid}.weka.trainingDataset.add( wekaInstance );
                            end
                            
                        end

                       waitbar(i/numel(annotationData.cellStats), hWaitBar, ...
                                  sprintf('Computing Features for file %d/%d in root directory %d/%d', ...
                                           fid, numel(annotationFileList), rid, numel(rootDirList) ));
                                   
                   end  

                   numWellSegmentedCells = numel(find(flagIsCellWellSegmented))
                   numCellsOnBorder = numel(find(flagIsCellOnBorder & flagIsCellWellSegmented))

                   for tid = 1:numel(classificationTask)

                       % write the features to disk
                       WriteFeatureMatrixToCSVFile( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s_All.csv', classificationTask{tid})), ...
                                                    cellFeatureData{tid}.featureMatrix, ...
                                                    cellFeatureData{tid}.featureLabels );

                       WriteFeatureMatrixToCSVFile( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s.csv', classificationTask{tid})), ...
                                                    cellFeatureData{tid}.featureMatrix(:, cellFeatureData{tid}.flagUseAsFeature), ...
                                                    cellFeatureData{tid}.featureLabels(cellFeatureData{tid}.flagUseAsFeature) );

                       annotationFilePath = curAnnotationFilePath;
                       featureData = cellFeatureData{tid};
                       featureData.featureFunctionName = sprintf( 'ComputeCellPatternFeatures_%s.m', ...
                                                                  classificationTask{tid} );                                                                  
                       save( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s.mat', classificationTask{tid})), ...
                             'annotationFilePath', ...
                             'featureData' );
                       
                   end

                   timeElapsed = toc

                end

            end       
            
            % write global files
            for tid = 1:numel(classificationTask)

               fprintf('\n\nClass Distribution for classification task - %s\n\n', classificationTask{tid} );

               classcolid = find(strcmpi(globalFeatureData{tid}.featureLabels, 'label.className'));
               tabulate( globalFeatureData{tid}.featureMatrix(:, classcolid) )

               WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, sprintf('cellCycleFeatures_%s_All.csv', classificationTask{tid})), ...
                                            globalFeatureData{tid}.featureMatrix, ...
                                            globalFeatureData{tid}.featureLabels );

               WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, sprintf('cellCycleFeatures_%s.csv', classificationTask{tid})), ...
                                            globalFeatureData{tid}.featureMatrix(:, globalFeatureData{tid}.flagUseAsFeature), ...
                                            globalFeatureData{tid}.featureLabels(globalFeatureData{tid}.flagUseAsFeature) );
               
               arffFileSaver = ArffSaver();
               arffFileName = java.lang.String( fullfile(outputRootDir, sprintf('cellCycleFeatures_%s.arff', classificationTask{tid})) );
               arffFileSaver.setFile( java.io.File(arffFileName) );
               arffFileSaver.setInstances( globalFeatureData{tid}.weka.trainingDataset );
               arffFileSaver.writeBatch();
                           
               annotationRootDirList = rootDirList;
               featureData = globalFeatureData{tid};
               featureData.featureFunctionName = sprintf( 'ComputeCellPatternFeatures_%s.m', ...
                                                          classificationTask{tid} );                     
               save( fullfile(outputRootDir, sprintf('cellCycleFeatures_%s.mat', classificationTask{tid})), ...
                     'annotationRootDirList', ...
                     'featureData' ); 
                 
            end

            % build classifiers for each stage 
            numEnsembleComponents = 3;
            
                % build a classifier for G1_S_G2_M     
                fprintf('\n\nBuilding G1_S_G2_M classifier ...\n\n' );
                
                tic                
                G1_S_G2M_Classifier = Vote();
                G1_S_G2M_Classifier.setCombinationRule( SelectedTag(Vote.MAJORITY_VOTING_RULE, Vote.TAGS_RULES) );
                
                for i = 1:numEnsembleComponents
                    
                    % create spread subsample filter
                    classBalancingFilter = SpreadSubsample();                    
                    classBalancingFilter.setDistributionSpread( 1.0 );
                    classBalancingFilter.setRandomSeed( i );
                    
                    % create libsvm classifier
                    svmClassifier = javaObject('weka.classifiers.functions.LibSVM');
                    svmClassifier.setCost( 512.0 );
                    svmClassifier.setGamma( 0.125 );
                    
                    % create filtered classfier
                    balancedClassifier = FilteredClassifier();
                    balancedClassifier.setFilter( classBalancingFilter );
                    balancedClassifier.setClassifier( svmClassifier );
                    
                    balancedClassifier.buildClassifier( globalFeatureData{1}.weka.trainingDataset );
                    
                    % add to ensemble
                    G1_S_G2M_Classifier.addPreBuiltClassifier( balancedClassifier );
                    
                end
                
                G1_S_G2M_Classifier.buildClassifier( globalFeatureData{1}.weka.trainingDataset );
                
                header = Instances(globalFeatureData{1}.weka.trainingDataset, 0);
                SerializationHelper.writeAll( fullfile(outputRootDir, 'G1_S_G2M.model'), ...
                                              javaArray('java.lang.Object',  G1_S_G2M_Classifier, header));         
                timeElapsed = toc
                
                % build a classifier for G2_M
                fprintf('\n\nBuilding G2_M classifier ...\n\n' );
                
                tic
                G2_M_Classifier = Vote();
                G2_M_Classifier.setCombinationRule( SelectedTag(Vote.MAJORITY_VOTING_RULE, Vote.TAGS_RULES) );
                
                for i = 1:numEnsembleComponents
                    
                    % create spread subsample filter
                    classBalancingFilter = SpreadSubsample();                    
                    classBalancingFilter.setDistributionSpread( 1.0 );
                    classBalancingFilter.setRandomSeed( i );
                    
                    % create libsvm classifier
                    randomForestClassifier = RandomForest();
                    
                    % create filtered classfier
                    balancedClassifier = FilteredClassifier();
                    balancedClassifier.setFilter( classBalancingFilter );
                    balancedClassifier.setClassifier( randomForestClassifier );
                    
                    balancedClassifier.buildClassifier( globalFeatureData{2}.weka.trainingDataset );
                    
                    % add to ensemble
                    G2_M_Classifier.addPreBuiltClassifier( balancedClassifier );
                    
                end
                
                G2_M_Classifier.buildClassifier( globalFeatureData{2}.weka.trainingDataset ); 
                
                header = Instances(globalFeatureData{2}.weka.trainingDataset, 0);
                SerializationHelper.writeAll( fullfile(outputRootDir, 'G2_M.model'), ...
                                              javaArray('java.lang.Object', G2_M_Classifier, header) );         
                timeElapsed = toc
                
            % cross validate models
                
                % G1_S_G2M classifier cross-validation
                fprintf('\n\nCross-validating G1_S_G2M classifier ...\n\n' );
                
                tic
                evaluator = Evaluation();
                evaluator.crossValidateModel(G1_S_G2M_Classifier, ...
                                             globalFeatureData{1}.weka.trainingDataset, ...
                                             Random(1));
    
                fprintf('\nResults Summary\n\n%s\n', char(evaluator.toSummaryString()) );
                fprintf('\nDetailed Accuracy By Class\n\n%s\n', char(evaluator.toClassDetailsString()) );
                fprintf('\nConfusion Matrix\n\n%s\n', char(evaluator.toMatrixString()) );
                
                GMean = 1;
                numClasses = globalFeatureData{1}.weka.trainingDataset.numClasses();
                for cid = 1:numClasses
                    GMean = GMean * evaluator.recall(cid);
                end
                GMean = GMean^(1/numClasses);
                fprintf('\nG-Mean - %f\n', G-Mean );                
                timeElapsed = toc
                
                % G1_S_G2M classifier cross-validation
                fprintf('\n\nCross-validating G1_S_G2M classifier ...\n\n' );
                
                tic
                evaluator = Evaluation();
                evaluator.crossValidateModel(G2_M_Classifier, ...
                                             globalFeatureData{2}.weka.trainingDataset, ...
                                             Random(1));
    
                fprintf('\nResults Summary\n\n%s\n', char(evaluator.toSummaryString()) );
                fprintf('\nDetailed Accuracy By Class\n\n%s\n', char(evaluator.toClassDetailsString()) );
                fprintf('\nConfusion Matrix\n\n%s\n', char(evaluator.toMatrixString()) );
                
                GMean = 1;
                numClasses = globalFeatureData{2}.weka.trainingDataset.numClasses();
                for cid = 1:numClasses
                    GMean = GMean * evaluator.recall(cid);
                end
                GMean = GMean^(1/numClasses);
                fprintf('\nG-Mean - %f\n', G-Mean );                
                timeElapsed = toc
                
            % store the models
            obj.G1_S_G2M_Model.classifier = G1_S_G2M_Classifier;
            obj.G1_S_G2M_Model.header = Instances(globalFeatureData{1}.weka.trainingDataset, 0);

            obj.G2_M_Model.classifier = G2_M_Classifier;
            obj.G2_M_Model.header = Instances(globalFeatureData{2}.weka.trainingDataset, 0);
            obj.flagModelReady = true;
            
            % switch off diary
            diary off;
            
        end % BuildModel
        
        function [ predictedLabels ] = predict(obj, imageData, imLabelCellSeg, cellStats)
        
            predictedLabels = cell( numel(cellStats), 1 );
            
            % weka imports
            import weka.core.*;
            import weka.core.converters.*;
            import weka.classifiers.*;
            
           % standardize the image data before computing any features
           imageDataAdjusted = StandardizeIntravitalImageData(imageData);

           % denoise it a bit using a median filter
           for i = 1:numel(imageData)
               imageDataAdjusted{i} = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted{i} );    
           end

           % compute features for each cell and run classifiers
           h = waitbar(0, 'Running Cell State Classifier ...' );
           
           flagStage2Tested = false;
           
           for i = 1:numel(cellStats)

                curCellStats = cellStats(i);

                fprintf( '\n>> Computing Features for cell %d/%d ...\n', i, numel(cellStats) );       

                % note down some cell indentity information 
                curCellStats.cellId = i;        
                imCurCellSegMask = (imLabelCellSeg == i);

                % Compute features for G1_S_G2M and run classifier
                [classNameList, featureStruct, ...
                featureVec, featureLabels, ...
                flagLearnFeature] = ComputeCellPatternFeatures_G1_S_G2M( imageDataAdjusted, ...
                                                                         imCurCellSegMask, ...
                                                                         curCellStats, false );
                
                [wekaInstance, ...
                 wekaAttributeList, ...
                 classIndex] = CellCycleClassifier.ConstructWekaInstance( featureVec(flagLearnFeature > 0), ...
                                                                          featureLabels(flagLearnFeature > 0), ...
                                                                          classNameList);
                                                                
                curDataSet = Instances('TestInstance', wekaAttributeList, 1 );
                curDataSet.add(wekaInstance);
                curDataSet.setClassIndex(classIndex-1);
                
                if i == 1 && ~obj.G1_S_G2M_Model.header.equalHeaders(curDataSet)
                    error( 'ERROR: mismatch between computed features and classifier features' );
                end
                
                predictedLabels{i} = classNameList{obj.G1_S_G2M_Model.classifier.classifyInstance(curDataSet.instance(0)) + 1};
                                                  
                % Compute features for G2_M
                if strncmpi( predictedLabels{i}, 'G2', 2 )
                    
                    [classNameList, featureStruct, ...
                     featureVec, featureLabels, ...
                     flagLearnFeature] = ComputeCellPatternFeatures_G2_M( imageDataAdjusted, ...
                                                                          imCurCellSegMask, ...
                                                                          curCellStats, false );

                    [wekaInstance, ...
                     wekaAttributeList, ...
                     classIndex] = CellCycleClassifier.ConstructWekaInstance( featureVec(flagLearnFeature > 0), ...
                                                                              featureLabels(flagLearnFeature > 0), ...
                                                                              classNameList);
                                                                          
                    curDataSet = Instances('TestInstance', wekaAttributeList, 1 );
                    curDataSet.add(wekaInstance);
                    curDataSet.setClassIndex(classIndex-1);
                    
                    if ~flagStage2Tested 
                        if ~obj.G2_M_Model.header.equalHeaders(curDataSet)                        
                            error( 'ERROR: mismatch between computed features and classifier features' );
                        else
                            flagStage2Tested = true;
                        end
                    end
                    
                    predictedLabels{i} = classNameList{obj.G2_M_Model.classifier.classifyInstance(curDataSet.instance(0)) + 1};
                    
                end

                waitbar(i/numel(cellStats), h, 'Running Cell State Classifier ...' );
           end              
            
           closeStatusDialog(h);
           
        end % predict
        
    end % methods
    
    methods (Static = true, Access = private)
       
        function [wekaInstance, wekaAttributeList, classIndex] = ConstructWekaInstance(featureVec, featureLabels, classNameList)

            wekaClassNameList = weka.core.FastVector();
            for i = 1:numel(classNameList)
                wekaClassNameList.addElement( java.lang.String(classNameList{i}) );
            end
            
            wekaAttributeList = weka.core.FastVector();                                
            wekaFeatureVec = zeros(1, numel(featureLabels));
            
            for attid = 1:numel(featureLabels)                                                                        

                curAttributeName = java.lang.String( featureLabels{attid} );

                if strcmpi(curAttributeName, 'label.className')
                    classIndex = attid;
                    curAttribute = weka.core.Attribute(curAttributeName, wekaClassNameList);
                    curAttributeVal = curAttribute.indexOfValue(java.lang.String(featureVec{attid}));                                        
                else                                     
                    curAttribute = weka.core.Attribute(curAttributeName);
                    curAttributeVal = featureVec{attid};
                end

                wekaAttributeList.addElement(curAttribute);
                wekaFeatureVec(attid) = curAttributeVal;

            end

            wekaInstance = weka.core.DenseInstance(1.0, wekaFeatureVec);
            
        end
        
    end
            
end % class def