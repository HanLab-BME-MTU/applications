
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    classificationTask = { 'G1_S_G2M', 'G2_M' };

    defaultDataRootDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis\annotations\cell_pattern_annotations';    
    defaultOutputDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis';    
    
    flagWriteGlobalFeatureFile = true;
    flagIgnoreCellsOnBorder = true;
    flagGenerateMergedFeatureFiles = true;
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ask the user to provide the list of directories to be processed
rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

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
    annotationFileList = rdir( fullfile(dataRootDir, '**', 'CellPatternAnnotation.mat') );

    % generate features for each annotation file
    globalFeatureData = cell(size(classificationTask));

    for fid = 1:numel(annotationFileList)

       fprintf( '\n\nComputing Features for file %d/%d ...\n\n', fid, numel(annotationFileList) );       

       tic

       curAnnotationFilePath = annotationFileList(fid).name   
       curAnnotationRelativeFilePath = annotationFileList(fid).name((numel(dataRootDir)+2):end);

       [pathstr, name, ext] = fileparts(curAnnotationRelativeFilePath); 
       curOutputDir = fullfile( outputDir, pathstr );
       if ~isdir(curOutputDir)
           mkdir(curOutputDir);
       end

       % load annotation data
       try
            annotationData = load( curAnnotationFilePath ); 
       catch
           fprintf('\nERROR: could not load annotation file\n' );
           continue;
       end

       % standardize the image data before computing any features
       imageDataAdjusted = cell(size(annotationData.imageData));
       for i = 1:numel(annotationData.imageData)
           imageDataAdjusted{i} = mat2gray(annotationData.imageData{i}) * 4096;
       end
       
       % denoise it a bit using a median filter
       parfor i = 1:numel(annotationData.imageData)
           imageDataAdjusted{i} = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted{i} );    
       end

       % make a note of all the cells touching the border
       borderCellIds = setdiff( unique(annotationData.imLabelCellSeg .* ~annotationData.imRegValidMask), 0 );

       % compute features for each cell
       cellFeatureData = cell(size(classificationTask));
       
       flagIsCellWellSegmented = zeros(numel(annotationData.cellStats),1);
       flagIsCellOnBorder = zeros(numel(annotationData.cellStats),1);
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

            end

       end  

       numWellSegmentedCells = numel(find(flagIsCellWellSegmented))
       numCellsOnBorder = numel(find(flagIsCellOnBorder & flagIsCellWellSegmented))

       for tid = 1:numel(classificationTask)
           if flagWriteGlobalFeatureFile
                if isempty(globalFeatureData{tid})
                    globalFeatureData{tid} = cellFeatureData{tid};
                else
                    globalFeatureData{tid}.featureStruct = cat(1, globalFeatureData{tid}.featureStruct, cellFeatureData{tid}.featureStruct);
                    globalFeatureData{tid}.featureMatrix = cat(1, globalFeatureData{tid}.featureMatrix, cellFeatureData{tid}.featureMatrix);
                end
           end

           % write the features to disk
           WriteFeatureMatrixToCSVFile( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s_All.csv', classificationTask{tid})), ...
                                        cellFeatureData{tid}.featureMatrix, ...
                                        cellFeatureData{tid}.featureLabels );

           WriteFeatureMatrixToCSVFile( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s.csv', classificationTask{tid})), ...
                                        cellFeatureData{tid}.featureMatrix(:, cellFeatureData{tid}.flagUseAsFeature), ...
                                        cellFeatureData{tid}.featureLabels(cellFeatureData{tid}.flagUseAsFeature) );

           annotationFilePath = curAnnotationFilePath;
           featureData = cellFeatureData{tid};
           save( fullfile(curOutputDir, sprintf('cellCycleFeatures_%s.mat', classificationTask{tid})), ...
                 'annotationFilePath', ...
                 'featureData' );
       end

       timeElapsed = toc

    end

    if flagWriteGlobalFeatureFile

       for tid = 1:numel(classificationTask)

           fprintf('\n\nClass Distribution for classification task - %s\n\n', classificationTask{tid} );

           classcolid = find(strcmpi(globalFeatureData{tid}.featureLabels, 'label.className'));
           tabulate( globalFeatureData{tid}.featureMatrix(:, classcolid) )

           WriteFeatureMatrixToCSVFile( fullfile(outputDir, sprintf('cellCycleFeatures_%s_All.csv', classificationTask{tid})), ...
                                        globalFeatureData{tid}.featureMatrix, ...
                                        globalFeatureData{tid}.featureLabels );

           WriteFeatureMatrixToCSVFile( fullfile(outputDir, sprintf('cellCycleFeatures_%s.csv', classificationTask{tid})), ...
                                        globalFeatureData{tid}.featureMatrix(:, globalFeatureData{tid}.flagUseAsFeature), ...
                                        globalFeatureData{tid}.featureLabels(globalFeatureData{tid}.flagUseAsFeature) );

           annotationFilePath = curAnnotationFilePath;
           featureData = globalFeatureData{tid};
           save( fullfile(outputDir, sprintf('cellCycleFeatures_%s.mat', classificationTask{tid})), ...
                 'annotationFilePath', ...
                 'featureData' );
       end

    end

    % switch off diary
    diary off;
    
end


% merge csv files pairwise
if flagWriteGlobalFeatureFile && flagGenerateMergedFeatureFiles

    filecombs = nchoosek( rootDirList, 2 );

    for i = 1:size(filecombs,1)

        fprintf('\n\nGenerating Merged Feature File - %d/%d\n\n', i, size(filecombs,1));

        [pathstr, dirName1, ~] = fileparts(filecombs{i,1});
        [pathstr, dirName2, ~] = fileparts(filecombs{i,2});

        for cid = 1:numel(classificationTask)

            % copy individual csv files to output directory
            featureFile1_src = fullfile(outputRootDir, dirName1, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
            featureFile1_dest = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s.csv', classificationTask{cid}, dirName1));
            copyfile(featureFile1_src, featureFile1_dest, 'f');
            ConvertCSVFileToArffFile( featureFile1_dest );

            featureFile2_src = fullfile(outputRootDir, dirName2, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
            featureFile2_dest = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s.csv', classificationTask{cid}, dirName2));
            copyfile(featureFile2_src, featureFile2_dest, 'f');      
            ConvertCSVFileToArffFile( featureFile2_dest );

            % append the two files and copy to output directory
            outFile = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s_%s.csv', classificationTask{cid}, dirName1, dirName2));
            AppendTwoFeatureFiles( featureFile1_src, featureFile2_src, outFile );
            ConvertCSVFileToArffFile( outFile );

        end

    end    

end
    
