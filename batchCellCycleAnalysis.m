
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    projectRootDir = '/home/drc16/intravital/CellCycleAnalysis';
    projectRootDirScratch = '/hms/scratch1/drc16/intravital/CellCycleAnalysis';
    
    % specify where the image data stored?
    dataRootDir = fullfile(projectRootDir, 'data' );            
    
    % specify where do you want the results to be stored?
    resultsRootDir = fullfile(projectRootDirScratch, 'results', 'batchAnalysis_322_M15');
    
    % specify the path to the region merging model file
    regionMergingModelFile = fullfile( projectRootDir, 'models', 'regionMerging', ...
                                       'train_M04_M09_M12', ...
                                       'regionMerging.model' );                                      
                                   
    % specify the path to the cell cycle state idenfitification model
    cellCycleModelFile = fullfile( projectRootDir, 'models', 'cellCycleClassification', ...
                                   'train_M04_M09_M12', ...
                                   'G1_S_G2_M.model' );
    
    % specify the path to the inventory file created in microsoft excel
    % each row of this file corresponds to one dataset 
    % it must contain columns containing: 
    % - path to file (*.oif, *.oib) containing the nuclear marker data
    % - path to file (*.oif, *.oib) containing the fucci cell cycle reporter
    % - mouse ID
    % - grid position indicating the locating of the tumor image
    %
    % in addition the inventory file can contain any columns that you want
    % to be passed on to the result file outputed
    inventoryFile = fullfile( resultsRootDir, 'analysis_file_for_batch_binned.xlsx' );
    
    nuclearMarkerFileColumnName = 'nuclear file'; 
    fucciFileColumnName = 'fucci file';
    mouseIDColumnName = 'Mouse';
    gridLocationColumnName = 'grid position';
    dataRootPrefix = 'Z:\intravital\data'; % to calculcate the relative path to result dir of each dataset
    defaultVoxelSpacing = [0.497, 0.497, 2]; % use this if metadata is missing or corrupted
    
    % run mode
    % testSingle - test run the analysis on a single dataset
    % testDeployAsJobs - test deploying the analysis a small set of data the cluster
    % deployAsJobs - deploy the analysis onto the cluster
    % collect - assembles per-dataset result files into one global file
    runModeOptions = {'testSingle', 'testDeployAsJobs', 'deployAsJobs', 'collect'};
    runMode = 3;
    
    % specify cluster resource requirements and submission queue
    numCoresRequested = [];
    memoryAmountInMB = 2^13;
    strSubmissionQueue = 'danuser_12h';
    
    % specify if you want to process specific files that failed earlier
    % empty recommended. Should have a very good reason to specify this
    fidRequestedAnalysisList = []; 
    
    flagCollectResultFilesAfterAnalysis = true;
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read and parse inventory file
PrettyPrintStepDescription( 'Reading Inventory file' );

inventoryFile

[numData, txtData, rawData] = xlsread(inventoryFile);
inventoryFileHeader = (rawData(1,:));

inventoryFileHeader'

% get the list of histone and fucci files and other columns needed
nuclearMarkerColumnId = find( strcmpi(inventoryFileHeader, nuclearMarkerFileColumnName) )
rawData(2:end, nuclearMarkerColumnId) = strtrim( strrep( strrep(rawData(2:end, nuclearMarkerColumnId), dataRootPrefix, ''), '\', '/') );

fucciColumnId = find( strcmpi(inventoryFileHeader, fucciFileColumnName) )
rawData(2:end, fucciColumnId) = strtrim( strrep( strrep(rawData(2:end, fucciColumnId), dataRootPrefix, ''), '\', '/' ) );

mouseColumnId = find( strcmpi(inventoryFileHeader, mouseIDColumnName) )
gridLocationColumnId = find( strcmpi(inventoryFileHeader, gridLocationColumnName) )
rawData(2:end, gridLocationColumnId) = strtrim( rawData(2:end, gridLocationColumnId) );

metaDataColumnIdList = setdiff(1:size(rawData,2), [nuclearMarkerColumnId, fucciColumnId]);
metaDataHeader = inventoryFileHeader( metaDataColumnIdList );

rawData = rawData(2:end,:);
numFiles = size(rawData,1)

% run modes
strRunMode = runModeOptions{runMode}

switch strRunMode

    case 'testSingle'

        PrettyPrintStepDescription( 'Testing batch analysis script' );
        
        fid = 1;
        curResultsRootDir = fullfile(resultsRootDir, 'testSingle');
        
        curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
        [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );
        
        curFucciFilePath = fullfile(dataRootDir, rawData{fid, fucciColumnId});

        curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
        curStrGridLocation = rawData{fid, gridLocationColumnId};
        curNuclearMarkerFileName
        curOutDir = fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName)        
        
        curMetaInfoStruct.header = metaDataHeader;
        curMetaInfoStruct.data = rawData(fid, metaDataColumnIdList);
        
        flagSuccess = performCellCycleAnalysis( curNuclearMarkerFilePath, curFucciFilePath, ...
                                                regionMergingModelFile, cellCycleModelFile, curOutDir, ...
                                                'defaultVoxelSpacing', defaultVoxelSpacing, ...
                                                'metaInfoStruct', curMetaInfoStruct, ... 
                                                'flagSaveImages', true );

        if ~flagSuccess
            error( 'Analysis Failed' );
        end
        
    case { 'testDeployAsJobs', 'deployAsJobs'}        

        PrettyPrintStepDescription( 'Deploying jobs onto the cluster' );
        
        if strcmp(runModeOptions{runMode}, 'testDeployAsJobs')
            numFiles = 3;
            curResultsRootDir = fullfile(resultsRootDir, 'testDeploy');            
        else
            curResultsRootDir = fullfile(resultsRootDir, 'analysis');
        end

        if ~isdir( curResultsRootDir )
            mkdir( curResultsRootDir);
        else
            if ~strcmp(input('Specified Result Directory already exists. Do you want to overwrite (y/n)?', 's'), 'y')
                return;                
            end
        end        
        
        % prepare a list of file IDs to be analyzed
        if strcmp(runModeOptions{runMode}, 'testDeployAsJobs') || isempty( fidRequestedAnalysisList )
            fidAnalysisList = 1:numFiles;
        else
            if any( fidRequestedAnalysisList < 1 | fidRequestedAnalysisList > numFiles )
                error( 'file ids must be in the range %d-%d', 1, numFiles );
            end
            fidAnalysisList = fidRequestedAnalysisList;
        end

        % setup job scheduler
        jm = findResource('scheduler','type','lsf');
        set(jm, 'ClusterMatlabRoot','/opt/matlab');
        
        strSubmitArguments = sprintf('-R "rusage[mem=%d]" -q %s', ...
                                      memoryAmountInMB, ...
                                      strSubmissionQueue )
                                 
        set(jm, 'SubmitArguments', strSubmitArguments );

        % setup job log data dir
        jobdatadir = fullfile(curResultsRootDir, 'jobdata');
        if isdir(jobdatadir)
            rmdir(jobdatadir, 's');
        end
        mkdir(jobdatadir);        
        set(jm, 'DataLocation', jobdatadir );

        % add a task for processing each file
        fprintf( '\n>> Creating and submitting jobs for processing each of the %d files ... \n', numel(fidAnalysisList) );
        
        taskList = cell(numFiles, 1);
        jobList = cell(numFiles, 1);
        
        statusOutDir = fullfile(curResultsRootDir, 'jobProgress');
        if isdir(statusOutDir)            
            rmdir(statusOutDir, 's');
        end
        mkdir(statusOutDir);
        
        for fid = fidAnalysisList

            curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
            [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );

            curFucciFilePath = fullfile(dataRootDir, rawData{fid, fucciColumnId});

            curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
            curStrGridLocation = rawData{fid, gridLocationColumnId};

            fprintf( '\n%.3d/%.3d: M%.2d -- %s -- %s \n', fid, numFiles, ...
                     rawData{fid, mouseColumnId}, curStrGridLocation, curNuclearMarkerFileName );            
            
            curOutDir = fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName);

            curMetaInfoStruct.header = metaDataHeader;
            curMetaInfoStruct.data = rawData(fid, metaDataColumnIdList);

            curFinishStatusReportFile = fullfile(statusOutDir, num2str(fid));
            
            funcArgs = { curNuclearMarkerFilePath, curFucciFilePath, ...
                         regionMergingModelFile, cellCycleModelFile, curOutDir, ...
                         'defaultVoxelSpacing', defaultVoxelSpacing, ...
                         'metaInfoStruct', curMetaInfoStruct, ... 
                         'finishStatusReportFile', curFinishStatusReportFile, ...
                         'flagSaveImages', true, 'flagParallelize', false };

            jobList{fid} = createJob(jm);            
            taskList{fid} = createTask(jobList{fid}, @performCellCycleAnalysis, 1, funcArgs, 'CaptureCommandWindowOutput', true );
	    submit(jobList{fid});
            
        end

        % submit job to cluster 
        fprintf( '\nAll of the %d jobs have been submitted successfully ...\n', numel(fidAnalysisList) );
        fprintf( '\n>>Now please wait for the job to finish ...\n' );
        
        jobTimer = tic;

        for i = 1:numel(fidAnalysisList)
            fid = fidAnalysisList(i);
            waitForState(jobList{fid});
            fprintf('\nFinished processing file %d/%d (%.2f%%)\n', fid, numFiles, 100.0 * i / numel(fidAnalysisList));  
        end
        
        totalJobTime = toc(jobTimer);        
        fprintf( '\nIt took a total of %f minutes to process %d files \n', totalJobTime / 60.0, numel(fidAnalysisList));
        
        % check if everything went well
        PrettyPrintStepDescription( 'Checking if everything went well' );

        flagTaskSuccess = false(numFiles,1);
        
        if isempty(fidRequestedAnalysisList)
            fidFail = fopen( fullfile(curResultsRootDir, 'analysisReport.log'), 'w' );
        else
            fidFail = fopen( fullfile(curResultsRootDir, 'analysisReportSpecific.log'), 'w' );
            fprintf(fidFail, '\nRequested analysis file id list: [%s]\n', ...
                             sprintf(' %d ', fidAnalysisList) );            
        end
        
        fprintf(fidFail, '\nIt took a total of %f minutes to process %d files \n', toc(jobTimer) / 60.0, numel(fidAnalysisList));

        for fid = fidAnalysisList

            curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
            [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );

            curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
            curStrGridLocation = rawData{fid, gridLocationColumnId};
            curOutDir = fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName);

            curTaskOutput = get(taskList{fid}, 'OutputArguments');

            if ~strcmp(taskList{fid}.State, 'failed') && ~isempty(curTaskOutput)
                flagTaskSuccess(fid) = curTaskOutput{1};
            end

            if ~flagTaskSuccess(fid)

                fprintf( '\n>> Failed to process file %d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );
                fprintf( fidFail, '\n>> Failed to process file %d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );
                
                if isempty( taskList{fid}.pGetError() )
                    strErrorMessage = sprintf( 'Failed file due to unknown error.' );
                else
                    strErrorMessage = taskList{fid}.pGetErrorMessage();
                end
                
                strErrorMessage = sprintf( '\n%s\nStartTime: %s\nFinish Time: %s\n', ...
                                           strErrorMessage, ...   
                                           taskList{fid}.StartTime, ...
                                           taskList{fid}.FinishTime );

                fprintf( fidFail, '\n%s\n', strErrorMessage );
                
            end
            
            fidLogFile = fopen( fullfile(curOutDir, 'performCellCycleAnalysis.log'), 'w' );
            
            if fidLogFile > 0
                
                fprintf( fidLogFile, '%s', taskList{fid}.CommandWindowOutput );
                
                if ~flagTaskSuccess(fid)
                    
                    fprintf( fidFail, '\nFor a detailed error report check log file in the result folder:\n%s\n', curOutDir );

                    fprintf( fidLogFile, '\nError:\n\n%s\n', strErrorMessage );
                    if ~isempty( taskList{fid}.pGetError() )
                        fprintf( fidLogFile, '\nDetailed Error Report:\n\n%s\n', ...
                                             taskList{fid}.pGetError().getReport() );
                                         
                    end
                    
                end

                fclose(fidLogFile);
                
            end
            
        end

        if any( ~flagTaskSuccess(fidAnalysisList) )
        
            flagAnalysisTaskSuccess = flagTaskSuccess(fidAnalysisList);
            
            fprintf( '\n>> Failed to process a total of %d/%d (%.2f%%) files ... \n', ...
                     sum(~flagAnalysisTaskSuccess), numel(fidAnalysisList), ...
                     100.0 * mean(~flagAnalysisTaskSuccess) );            

            fprintf( fidFail, '\n>> Failed to process a total of %d/%d (%.2f%%) files ... \n', ...
                               sum(~flagAnalysisTaskSuccess), numel(fidAnalysisList), ...
                               100.0 * mean(~flagAnalysisTaskSuccess) );            
                 
            strFailedFileID = sprintf( ' %d ', fidAnalysisList((find(~flagAnalysisTaskSuccess))') );
            fprintf( '\n>> Failed file id list:\n\n [%s] \n', strFailedFileID );
            fprintf( fidFail, '\n>> Failed file id list:\n\n [%s] \n', strFailedFileID );
            
        else
            
            fprintf( '\nEurekaaaaa !!! All files were processed successfully ... \n' );
            fprintf( fidFail, '\nEurekaaaaa !!! All files were processed successfully ... \n' );
            
        end
        
        fclose(fidFail);
        
        % destroy the jobs
        for fid = fidAnalysisList
            destroy(jobList{fid});
        end
        
end

% collect result files
if strcmp(strRunMode, 'collect') || flagCollectResultFilesAfterAnalysis
    
    PrettyPrintStepDescription( 'Collecting result files' );

    curResultsRootDir = fullfile(resultsRootDir, 'analysis');
    globalStackAnalysisFile = fullfile(curResultsRootDir, 'stackAnalysisInfo.csv');
    globalCellAnalysisFile = fullfile(curResultsRootDir, 'cellAnalysisInfo.csv');
    fidResultCollection = fopen( fullfile(curResultsRootDir, 'resultFileCollection.log'), 'w' );
    
    fidAnalysisFileNotFound = [];
    
    for fid = 1:numFiles

        curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
        [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );

        fprintf( '\nCollecting analysis file %d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );

        curFucciFilePath = fullfile(dataRootDir, rawData{fid, fucciColumnId});

        curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
        curStrGridLocation = rawData{fid, gridLocationColumnId};

        curOutDir = fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName);

        curStackFile = fullfile(curOutDir, 'stackAnalysisInfo.csv');
        curCellFile = fullfile(curOutDir, 'cellAnalysisInfo.csv'); 

        if ~exist(curStackFile, 'file') || ~exist(curCellFile, 'file')
            
            fprintf( '\n\tUnable to find cell and stack analysis files\n' );
            fprintf( fidResultCollection, '\nAnalysis files not found for dataset %d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );
            fidAnalysisFileNotFound(end+1) = fid;
            continue;
        end
        
        if fid == 1
           copyfile( curStackFile, globalStackAnalysisFile, 'f');
           copyfile( curCellFile, globalCellAnalysisFile, 'f');
        else

            AppendTwoFeatureFiles(globalStackAnalysisFile, curStackFile);
            AppendTwoFeatureFiles(globalCellAnalysisFile, curCellFile);

        end

    end     
    
    if ~isempty(fidAnalysisFileNotFound)
        
        fprintf( '\nUnable to find cell and stack analysis for the following %d datasets:\n\n[%s]\n', ...
                 numel(fidAnalysisFileNotFound) , sprintf(' %d ', fidAnalysisFileNotFound) );

        fprintf( fidResultCollection, '\nUnable to find cell and stack analysis for the following %d datasets:\n\n[%s]\n', ...
                                       numel(fidAnalysisFileNotFound), ...
                                       sprintf(' %d ', fidAnalysisFileNotFound) );            
    else

        fprintf( fidResultCollection, '\nHurraayyy!!!\n\nValid result files were found and merged into the global file for all of the %d datsets\n', numFiles );            
        
    end
    
    fclose(fidResultCollection);
    
end
