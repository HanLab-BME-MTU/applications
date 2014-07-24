%function trainRegionMergingModel_3D

%     clc
%     clear 

    %***********************************************************************
    % 
    %                           PARAMETERS 
    %
    %***********************************************************************

        defaultDataRootDir = 'Z:\intravital\data';    
        defaultOutputDir = 'Z:\intravital\data';    

    %***********************************************************************

    % ask the user to provide the list of directories to be processed
    rootDirList = uipickfiles('Prompt', 'Select folders containing the annotation files', ...
                              'REFilter', 'DNADamageAnnotation\.mat$');

    if ~iscell(rootDirList)
        return;
    end
    
    % ask the user to select the output directory
    outputRootDir = uigetdir(defaultOutputDir, 'Select Output Directory' ); 

    if ~ischar(outputRootDir)
       return; 
    end
    
    hStatusDialog = waitbar(0, 'Building Model ...This might take a while, you may want to get some coffee.');
    
    % make a list of files that need to be processed
    annotationFileList = [];
    
    for rid = 1:numel(rootDirList)

        curAnnotationFileList = rdir( fullfile(rootDirList{rid}, '**', 'DNADamageAnnotation.mat') );
        annotationFileList = cat(1, annotationFileList, curAnnotationFileList );
        
    end
    
    if isempty(annotationFileList)        
        h = msgbox('No annotation files were found in the specified folders', ...
                   'Region merging model training', 'help', 'modal' );
        closeStatusDialog(hStatusDialog);
        return;
    end
    
    % log
    diary off;
    diary_file = fullfile( outputRootDir , sprintf( '%s.log' , mfilename ) );
    if exist( diary_file, 'file' )
        fclose all;
        delete( diary_file );
    end
    diary( diary_file );
    diary on;

    outputRootDir
    
    % process each annotation file
    numMerge = 0;
    numDontMerge = 0;
    
    featureData = [];
    featureClass = [];
        
    for fid = 1:numel(annotationFileList)

        fprintf( '\n\n*****************************************************\n\n' );           

        fprintf( '\n\nComputing Features for file %d/%d ...\n\n', fid, numel(annotationFileList) );       

        feaCompTimer = tic;

        curAnnotationFilePath = annotationFileList(fid).name   

        % load annotation data
        try
            annotationData = load( curAnnotationFilePath ); 
        catch
           fprintf('\nERROR: could not load annotation file\n' );
           continue;
        end

       % get data
       imInput = annotationData.imageData{annotationData.channelId53BP1};
       spacing = annotationData.metadata.voxelSpacing;       
       
       imLabelFociSeg = annotationData.imLabelFociSeg;
       imFociSeedPoints = annotationData.imFociSeedPoints;

       fociStats = annotationData.fociStats;       
       fociStats([fociStats.isValid] < 0) = []; % drop unannotated foci
       
       % compute foci detection features
       featureDataLocal = [];
       featureClassLocal = [];
       
       for i = 1:numel(fociStats)
           
           if fociStats(i).isValid > 0
               featureClassLocal = cat(1, featureClassLocal, {'Good_Detection'});               
           else
               featureClassLocal = cat(1, featureClassLocal, {'Bad_Detection'});
           end
           
           [featureDataStruct] = ComputeFociDetectionFeatures(imInput, fociStats, i, 'spacing', spacing);          
           [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureDataStruct );        
           featureDataLocal = cat(1, featureDataLocal, featureVec);
           
       end
       
       fprintf( '\nclass distribution:\n' );
       tabulate( featureClassLocal ) 
       
       featureData = cat(1, featureData, featureDataLocal);       
       featureClass = cat(1, featureClass, featureClassLocal);
       
       computationTime = toc(feaCompTimer)
       
       waitbar(fid / numel(annotationFileList), hStatusDialog);
       
    end

   fprintf( '\n\n*****************************************************\n\n' );
   
    if isempty(featureClass)        
        msgbox('No valid annotation files were found in the specified folders', ...
               'Region merging model training', 'error', 'modal' );
        closeStatusDialog(hStatusDialog);
        return;
    end
   
   fprintf( '\nTotal class distribution: \n' );
   tabulate( featureClass );
        
   % write csv file -- without debug info
    fprintf( '\nWriting csv file without region info ... \n' );
   
   featureLabels = cat(2, 'ClassLabel', featureNameList);
   featureMatrix = cat(2, featureClass, featureData);
   WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, 'fociDetectionFeatures.csv'), featureMatrix, featureLabels );

   % write arff file for weka
   fprintf( '\nConverting csv file to arff file ... \n' );
   
   AddWekaClassesToPath();

   ConvertCSVFileToArffFile( fullfile(outputRootDir, 'fociDetectionFeatures.csv') );
   
   % generate classification model
   wekaTrainingDataset = getWekaDatasetFromCSVFile( fullfile(outputRootDir, 'fociDetectionFeatures.csv') );
   wekaTrainingDataset.setClassIndex(0);
   wekaModel = BuildWekaBalancedRandomForestEnsemble(wekaTrainingDataset, fullfile(outputRootDir, 'fociDetection.model'));
   
   %wrapup
   closeStatusDialog(hStatusDialog);
   diary off;
   
%end
    
    
    
    
    
    
    
    
    