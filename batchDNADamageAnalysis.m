
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    projectRootDir = '/home/drc16/intravital/DNADamageAnalysis';
    projectRootDirScratch = '/hms/scratch1/drc16/intravital/DNADamageAnalysis';

    projectRootDir = 'Z:\intravital\DNADamageAnalysis';
    projectRootDirScratch = 'Z:\intravital\DNADamageAnalysis';
    
    % specify where the image data stored?
    dataRootDir = fullfile(projectRootDir, 'data' );            
    
     % specify where do you want the results to be stored?
    resultsRootDir = fullfile(projectRootDirScratch, 'results', 'batchAnalysis_M14_to_M17');
    
    % specify the path to the region merging model file
    regionMergingModelFile = fullfile( projectRootDir, 'models', 'regionMerging', ...
                                       'train_M04_M09_M12', ...
                                       'regionMerging.model' );                                      
    
    % specify the path to the inventory file created in microsoft excel
    % each row of this file corresponds to one dataset 
    % it must contain columns containing: 
    % - realtive path to file (*.oif, *.oib) containing the image data
    % - mouse ID
    % - 53BP1 channel index
    %
    % In addition the inventory file can contain any columns that you want
    % to be passed on to the result files outputed by the analysis scripts
    inventoryFile = fullfile( resultsRootDir, 'inventory_file.xlsx' );

    colName_RelImageFilePath = 'RelImageFilePath';
    colName_MouseId = 'MouseId';
    colName_53BP1ChannelId = '53BP1ChannelId';
    dataRootPrefix = ''; % to calculcate the relative path to result dir of each dataset
    
    % run mode
    % testSingle - test run the analysis on a single dataset
    % testDeploy - test deploying the analysis a small set of data the cluster
    % deploy - deploy the analysis onto the cluster
    % collect - assembles per-dataset result files into one global file
    runModeOptions = {'testSingle', 'testDeploy', 'deploy', 'collect'};
    runMode = 1;
    
    % Do you want to save images?
    flagSaveImages = true;
    
    % specify cluster resource requirements and submission queue
    numCoresRequested = [];
    memoryAmountInMB = 8 * 2^10; % ~8000MB
    strSubmissionQueue = 'danuser_12h';
    
    % specify if you want to process specific files that failed earlier
    % empty recommended. Should have a very good reason to specify this.
    % I use this if a small number of nodes of the cluster failed for 
    % some external reason while running the requested batch job.
    fidRequestedAnalysisList = []; 
    
    flagCollectResultFilesAfterAnalysis = true;
    resultFileList = {'stackAnalysisInfo.csv', 'cellAnalysisInfo.csv', 'fociAnalysisInfo.csv'};
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read and parse inventory file
PrettyPrintStepDescription( 'Reading Inventory file' );

inventoryFile

[numData, txtData, rawData] = xlsread(inventoryFile);
inventoryFileHeader = (rawData(1,:));

inventoryFileHeader'

% get the list of image data files and other columns needed
colId_RelImageFilePath = find( strcmpi(inventoryFileHeader, colName_RelImageFilePath) )
rawData(2:end, colId_RelImageFilePath) = strtrim( strrep( strrep(rawData(2:end, colId_RelImageFilePath), dataRootPrefix, ''), '\', '/') );
colId_Mouse = find( strcmpi(inventoryFileHeader, colName_MouseId) )
colId_53BP1ChannelId = find( strcmpi(inventoryFileHeader, colName_53BP1ChannelId) )

metaDataColumnIdList = 1:size(rawData,2);
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
        
        curImageFilePath = fullfile(dataRootDir, rawData{fid, colId_RelImageFilePath});
        [pname, curImageFileName, fext] = fileparts( curImageFilePath );
        
        curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, colId_Mouse} );
        cur53BP1ChannelId = rawData{fid, colId_53BP1ChannelId};
        
        curImageFileName
        
        curOutDir = fullfile(curResultsRootDir, curStrMouseID, curImageFileName)        
        
        curMetaInfoStruct.header = metaDataHeader;
        curMetaInfoStruct.data = rawData(fid, metaDataColumnIdList);
        
        flagSuccess = performDNADamageAnalysis( curImageFilePath, cur53BP1ChannelId, ...
                                                regionMergingModelFile, curOutDir, ...
                                                'metaInfoStruct', curMetaInfoStruct);

        if ~flagSuccess
            error( 'Analysis Failed' );
        end
        
    case { 'testDeploy', 'deploy'}        

        PrettyPrintStepDescription( 'Deploying jobs onto the cluster' );
        
        if strcmp(runModeOptions{runMode}, 'testDeploy')
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
        if strcmp(runModeOptions{runMode}, 'testDeploy') || isempty( fidRequestedAnalysisList )
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

            curImageFilePath = fullfile(dataRootDir, rawData{fid, colId_RelImageFilePath});
            [pname, curImageFileName, fext] = fileparts( curImageFilePath );

            curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, colId_Mouse} );
            cur53BP1ChannelId = rawData{fid, colId_53BP1ChannelId};
            
            fprintf( '\n%.3d/%.3d: M%.2d -- %s \n', fid, numFiles, ...
                     rawData{fid, colId_Mouse}, curImageFileName );            
            
            curOutDir = fullfile(curResultsRootDir, curStrMouseID, curImageFileName);

            curMetaInfoStruct.header = metaDataHeader;
            curMetaInfoStruct.data = rawData(fid, metaDataColumnIdList);

            curFinishStatusReportFile = fullfile(statusOutDir, num2str(fid));
            
            funcArgs = { curImageFilePath, cur53BP1ChannelId, ...
                         regionMergingModelFile, curOutDir, ...
                         'metaInfoStruct', curMetaInfoStruct, ... 
                         'finishStatusReportFile', curFinishStatusReportFile, ...
                         'flagSaveImages', true, 'flagParallelize', false };

            jobList{fid} = createJob(jm);            
            taskList{fid} = createTask(jobList{fid}, @performDNADamageAnalysis, 1, funcArgs, 'CaptureCommandWindowOutput', true );
            
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

            curImageFilePath = fullfile(dataRootDir, rawData{fid, colId_RelImageFilePath});
            [pname, curImageFileName, fext] = fileparts( curImageFilePath );

            curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, colId_Mouse} );
            curOutDir = fullfile(curResultsRootDir, curStrMouseID, curImageFileName);

            curTaskOutput = get(taskList{fid}, 'OutputArguments');

            if ~strcmp(taskList{fid}.State, 'failed') && ~isempty(curTaskOutput)
                flagTaskSuccess(fid) = curTaskOutput{1};
            end

            if ~flagTaskSuccess(fid)

                fprintf( '\n>> Failed to process file %d/%d: %s\n', fid, numFiles, curImageFileName );
                fprintf( fidFail, '\n>> Failed to process file %d/%d: %s\n', fid, numFiles, curImageFileName );
                
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
            
            fidLogFile = fopen( fullfile(curOutDir, 'performDNADamageAnalysis.log'), 'w' );
            
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

    if strcmp(strRunMode, 'collect')
        curResultsRootDir = fullfile(resultsRootDir, 'analysis');
    end
    
    globalStackAnalysisFile = fullfile(curResultsRootDir, 'stackAnalysisInfo.csv');
    globalCellAnalysisFile = fullfile(curResultsRootDir, 'cellAnalysisInfo.csv');
    globalFociAnalysisFile = fullfile(curResultsRootDir, 'fociAnalysisInfo.csv');
    fidResultCollection = fopen( fullfile(curResultsRootDir, 'resultFileCollection.log'), 'w' );
    
    fidAnalysisFileNotFound = [];
    
    for fid = 1:numFiles

        curImageFilePath = fullfile(dataRootDir, rawData{fid, colId_RelImageFilePath});
        [pname, curImageFileName, fext] = fileparts( curImageFilePath );

        fprintf( '\nCollecting analysis file %d/%d: %s\n', fid, numFiles, curImageFileName );

        curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, colId_Mouse} );

        curOutDir = fullfile(curResultsRootDir, curStrMouseID, curImageFileName);

        % check if all result files are present
        flagAllResultFilesFound = true;
        for afid = 1:numel(resultFileList)
        
            curLocalAnalysisFile = fullfile(curOutDir, resultFileList{afid});
            if ~exist(curLocalAnalysisFile, 'file')
                fprintf('\n\tCould not find analysis file %s in the result output directory %s\n', ...
                         resultFileList{afid}, curOutDir);
                fprintf(fidResultCollection, '\nAnalysis file %s not found for dataset %d/%d: \n%s\n', ...
                                             resultFileList{afid}, fid, numFiles, curImageFileName );
                flagAllResultFilesFound = false;
                fidAnalysisFileNotFound(end+1) = fid;
                break;
            end
            
        end
        
        if ~flagAllResultFilesFound
            continue;
        end
        
        % append to global file
        for afid = 1:numel(resultFileList)
            
            curLocalAnalysisFile = fullfile(curOutDir, resultFileList{afid});
            curGlobalAnalysisFile = fullfile(curResultsRootDir, resultFileList{afid});
            
            if fid == 1
               copyfile(curLocalAnalysisFile, curGlobalAnalysisFile, 'f');
            else
                AppendTwoFeatureFiles(curGlobalAnalysisFile, curLocalAnalysisFile);
            end
            
        end
        
    end     
    
    if ~isempty(fidAnalysisFileNotFound)
        
        fprintf( '\nUnable to find all analysis result files for the following %d datasets:\n\n[%s]\n', ...
                 numel(fidAnalysisFileNotFound) , sprintf(' %d ', fidAnalysisFileNotFound) );

        fprintf( fidResultCollection, '\nUnable to find cell and stack analysis for the following %d datasets:\n\n[%s]\n', ...
                                       numel(fidAnalysisFileNotFound), ...
                                       sprintf(' %d ', fidAnalysisFileNotFound) );            
    else

        fprintf( fidResultCollection, '\nHurraayyy!!!\n\nValid result files were found and merged into the global file for all of the %d datsets\n', numFiles );            
        
    end
    
    fclose(fidResultCollection);
    
end

