
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    projectRootDir = 'Z:\intravital';
    projectRootDirScratch = 'c:\deepak';

    % where is the image data stored?
    dataRootDir = fullfile(projectRootDir, 'data' );            
    
    % whare do you want the results to be stored?
    resultsRootDir = fullfile(projectRootDirScratch, 'results', 'batchAnalysisTest');
    
    % path to the region merging model file
    regionMergingModelFile = fullfile( projectRootDir, 'data', 'image_analysis', 'models', ...
                                       'region_merging', 'models', 'train_M04_M09_M12', ...
                                       'regionMerging.model' );                                      

    % path to the cell cycle state identification model file
    cellCycleModelFile = fullfile( projectRootDir, 'data', 'image_analysis', 'models', ...
                                   'cell_state_classification', 'models', 'train_M04_M09_M12', ...
                                   'G1_S_G2_M.model' );
    
    % path to inventory file created in microsoft excel
    inventoryFile = fullfile( projectRootDirScratch, 'results', 'batchAnalysis', 'inventory_files', 'analysis_file_for_batch.xlsx' );
    nuclearMarkerFileColumnName = 'nuclear file';
    fucciFileColumnName = 'fucci file';
    mouseIDColumnName = 'Mouse';
    gridLocationColumnName = 'grid position';
    dataRootPrefix = 'Z:\intravital\data';
    
    % run mode
    runModeOptions = {'test', 'test_deploy', 'deploy', 'collect'};
    runMode = 1;
    
    numCoresRequested = [];
    memoryAmountInMB = 2^13;
    strSubmissionQueue = 'danuser_2h';
    
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

% run analysis
curResultsRootDir = fullfile(resultsRootDir, 'analysis');

flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
if ~flagPoolOpenedAlready 
    matlabpool open;
end            

flagTaskSuccess = false(numFiles, 1);

jobTimer = tic;

for fid = 1:numFiles

    PrettyPrintStepDescription( sprintf( 'Running cell cycle analysis for file %d/%d', fid, numFiles) );
    
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

    flagTaskSuccess(fid) = performCellCycleAnalysis( curNuclearMarkerFilePath, curFucciFilePath, ...
                                                     regionMergingModelFile, cellCycleModelFile, curOutDir, ...
                                                     'metaInfoStruct', curMetaInfoStruct, ... 
                                                     'flagSaveImages', true, 'flagParallelize', true  );
   
end

totalJobTime = toc(jobTimer);

% check for any errors
PrettyPrintStepDescription( 'Checking if everything went well' );

fprintf( '\nIt took a total of %f minutes to process %d files \n', totalJobTime / 60.0, numFiles );

if any( ~flagTaskSuccess )

    fidFail = fopen( fullfile(curResultsRootDir, 'analysisErrorInventory.log'), 'w' );

    fprintf(fidFail, '\nIt took a total of %f minutes to process %d files \n', totalJobTime / 60.0, numFiles );

    fprintf( '\n>> Failed to process %d/%d (%.2f%%) files ... \n', sum(~flagTaskSuccess), numFiles, 100.0 * mean(~flagTaskSuccess) );            
    fprintf( fidFail, '\n>> Failed to process %d/%d (%.2f%%) files ... \n', sum(~flagTaskSuccess), numFiles, 100.0 * mean(~flagTaskSuccess) );

    for fid = (find(~flagTaskSuccess))'

        curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
        [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );
        curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
        curStrGridLocation = rawData{fid, gridLocationColumnId};
        
        curOutDir = fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName);

        fidLogFile = fopen( fullfile(curOutDir, 'performCellCycleAnalysis.log'), 'a' );
        
        fprintf( '\n%d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );
        fprintf( fidFail, '\n>> %d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );

        fclose(fidLogFile);

    end

    fclose(fidFail);

else

    fprintf( '\nAll files were processed success fully ... \n' );

end

% now assemble the csv files
PrettyPrintStepDescription( 'Collecting result files' );

globalStackAnalysisFile = fullfile(resultsRootDir, 'analysis', 'stackAnalysisInfo.csv');
globalCellAnalysisFile = fullfile(resultsRootDir, 'analysis', 'cellAnalysisInfo.csv');

last_percent_done = 0;
numPrint = 0;
flagFirstFile = true;

for fid = (find(flagTaskSuccess))'

    curNuclearMarkerFilePath = fullfile(dataRootDir, rawData{fid, nuclearMarkerColumnId});
    [pname, curNuclearMarkerFileName, fext] = fileparts( curNuclearMarkerFilePath );

    curFucciFilePath = fullfile(dataRootDir, rawData{fid, fucciColumnId});

    curStrMouseID = sprintf('Mouse_%.2d', rawData{fid, mouseColumnId} );
    curStrGridLocation = rawData{fid, gridLocationColumnId};

    curOutDir = strtrim( fullfile(curResultsRootDir, curStrMouseID, curStrGridLocation, curNuclearMarkerFileName) );

    fprintf( '\n%d/%d: %s\n', fid, numFiles, curNuclearMarkerFileName );
    
    curStackFile = fullfile(curOutDir, 'stackAnalysisInfo.csv');
    curCellFile = fullfile(curOutDir, 'cellAnalysisInfo.csv'); 

    if flagFirstFile        
       copyfile( curStackFile, globalStackAnalysisFile, 'f');
       copyfile( curCellFile, globalCellAnalysisFile, 'f');
       flagFirstFile = false;
    else

        AppendTwoFeatureFiles(globalStackAnalysisFile, curStackFile);
        AppendTwoFeatureFiles(globalCellAnalysisFile, curCellFile);

    end

    percent_done = round(100*fid/numFiles);       

    if percent_done > last_percent_done
        fprintf( '%.2d%%  ', percent_done );
        last_percent_done = percent_done;
        numPrint = numPrint + 1;
        if mod( numPrint, 10 ) == 0
           fprintf( '\n' ); 
        end
    end        

end

