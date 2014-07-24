function [trajectoryDescription] = trajectoryAnalysis(inputData,ioOpt,testOpt)
%TRAJECTORYANALYSIS analyzes experimental and simulated microtubule trajectories and returns corresponding statistical descriptors
%
% SYNOPSIS  [trajectoryDescription] = trajectoryAnalysis(inputData,ioOpt,testOpt)
%
% INPUT  inputData(1:n) (opt): structure containing n different trajectories with fields
%           - distance   tx2 array [distance, sigmaDistance] in microns
%           - time       tx2 array [time, sigmaTime] in seconds
%           - timePoints tx1 array [timePoint#] 
%           - info       struct containing additional information about the inputData
%                  -tags    : cell containing the two strings designating
%                             the tags between which the distance is
%                             measured, e.g. {'spb1','cen1'}
%
%       alternatively, inputData can be a collection of several sets of data:
%         inputData(1:iRun) with fields
%           - data              : (as inputData(1:n) above)
%           - fileNameList{1:n} : a cell array of strings of the same
%                                  length as data containing fileNames/descriptions 
%                                  of the individual trajectories           
%
%        ioOpt: optional structure with the following optional fields
%           - details       : return detailed analysis? [{0}/1]
%           - convergence   : add convergence data? [{0}/1]
%           - verbose       : whether to show any graphs/write to
%                             commandline or not. (vector)
%                             Verbose = 1 writes to commandline
%                             Verbose = 2 shows trajectories
%                             Verbose = 3 shows distributions 
%                                           (currently not supported)
%                             Verbose = 4 shows clustering result
%                             Verbose = 5 shows time distribution
%                             Verbose = [1,2,3,4,5] shows all
%                             Verbose = [] or [0] does not show anything.
%                             Asking for savePaths is not affected by
%                             verbose settings
%           - loadData      : [{0}/1/n] if inputData is not empty: whether to
%                               ask for more data or not. if inputData is
%                               not a run, loadData is forced to 1 (if n) and the
%                               data is appended to the input
%           - saveTxt       : write statistics in text file? [0/{1}] 
%                               currently works for 1 run only
%           - saveTxtPath   : if saveTxt is 1, specify the path to save.
%                               if does not end with filesep, program
%                               assumes it to include the filename!
%                               otherwise the program asks for the save-path
%           - saveMat       : save output in mat file? [{0}/1]
%                               will save all the runs in one single file
%           - saveMatPath   : if saveMat is 1, specify the path to save.
%                               if does not end with filesep, program
%                               assumes it to include the filename!
%                               otherwise the program asks for the save-path
%           - fileNameList  : list of fileNames/descriptions of the input
%                               data. This only needs to be specified if
%                               the inputData is not given as run
%           - expOrSim      : switch that tells whether inputData is from
%                               experiment or simulation ['e'/'s']. def: 'x'
%           - calc2d        : [{0}/1/2] depending on whether the
%                               normal 3D-data or the maxProjection
%                               or the in-focus-slice data should be
%                               used. Works only if inputData is being loaded by
%                               the program. Otherwise, you have to specify
%                               this option in
%                               calculateTrajectoryFromIdlist
%           - clusterData   : [{0}/1/2/3] uses EM clustering to find clusters
%                               of speeds. 1 clusters only at the end, 2
%                               will additionaly cluster for the convergence
%                               statistics, and 3 will cluster all
%           - realTime      : [0/{1}] whether or not to use the real
%                               timestamp. If 0, the intervals will be
%                               rounded to entire seconds. Warning: this
%                               will NOT affect how frequencies are being
%                               calculated!
%                               
%        testOpt: optional structure with the following optional fields
%           - splitData     : [{0}/1] split inputData into two sets, returns two
%                               output arguments
%           - subset        : [fraction,n] {[1,1]} randomly takes n unique
%                              subsets of size fraction*nData. MaxN is 100
%                              (exclusive to splitData).
%           - debug         : [{0}/1] turns on debug features
%           - randomize     : [{0}/1] randomize the input list
%           - clustermax    : max # of clusters the EM algorithm looks for
%           - downsample    : [{1}/n] takes only every nth point, returns n
%                               times more trajctories
%
% OUTPUT trajectoryDescription   nRuns by 1 structure with fields
%           .individualStatistics(1:n)  statistics of individual trajectory
%               .statistics             statistics (fieldnames: see below)
%               .details (opt)         
%                    .dataListGroup     fits to the trajectory
%                    .dataListSeed      single-interval classification
%                    .distributions     list of antipoleward&poleward
%                                       speeds, growt&shrinkage times,
%                                       distances, pause times.
%                                       [value,number of repeats,sigma]
%           .overallStatistics          mean over all trajectories in run
%               OR
%           .convergenceStatistics(1:n) mean over first i trajectories     
%           .overallDistribution        catenation of individual
%                                       distributions
%           .info                       only in first element of trajDes:
%                                       info on statistics and lists
%                                       ## not implemented yet ##
%
%        columns of dataLists
%           1:startIdx, 2:endIdx, 3:state, 4:slope, 5:slopeSigma, 6:slopeSigmaAPR,
%           7:deltaT, 8:(deltaTSigma), 9:deltaD, 10:deltaDSigma, 11:startDistance
%
%        fields of statistics-structure 
%           SEM: standard error of the mean
%           STD: standard deviation of the sample (SEM*sqrt(n))
%
%       WARNING: FREQUENCIES ARE CALCULATED AS PER TIMEPOINT. TO GET
%       FREQUENCIES PER SECOND, DIVIDE BY THE SAMPLING RATE IN SECONDS
%     'ap2tpFreq__cat' ,        catastrophe frequency (m,sem,n; [-])
%     'tp2apFreq__res' ,        rescue frequency      (m,sem,n; [-])
%     'antipolewardSpeed' ,     growth speed          (m,sem,n; [um/min])
%     'polewardSpeed' ,         shrinkage speed       (m,sem,n; [um/min])
%     'distanceMean',           mean spb-cen distance (m,sem; [um])
%     'distanceStd'             std of distance       (m,sem; [um])
%           for individual statistics, there can be no sem
%     'minDistance',            global minimum distance (m,std; [um])
%     'minDistanceM5' ,         mean of 5 smallest distances (m,sem; [um])
%     'maxDistance' ,           global maximum distance (m,std; [um])  
%     'maxDistanceM5' ,         mean of 5 largest distances (m,sem; [um])
%     'pauseNumber',            number of pause events
%     'avgApDistance' ,         mean distance per growth event (m,sem; [um])
%     'avgTpDistance' ,         mean distance per shrinkage event (m,sem; [um])
%     'avgUndetDistance' ,      avg of absolute distance in undet. intervals (m,sem; [um])
%     'antipolewardTime' ,      total AP time [s]; % of total traj. time     
%     'polewardTime' ,          total TP time [s]; % of total traj. time
%     'pauseTime' ,             total pause time [s]; % of total traj. time 
%     'undeterminedTime' ,      total undet. time [s]; % of total traj. time 
%     'deletedTime' ,           total not analyzed time [s] - not counting
%                               deletion at the end of the trajectory
%     'nTimepoints',            number of total timepoints; avg per trajectory
%     'averageTime',            average timeLapse [s], std
%        
%
%c: 11/03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%===================================================================
%--------------DEFINE CONSTANTS ET AL-- (for test input, see below)
%===================================================================

%set defaults
showDetails = 1;
doConvergence = 0;
saveTxt = 1;
saveTxtPath = '';
saveTxtName = '';
saveMat = 0;
saveMatPath = '';
saveMatName = '';
verbose = [1,2]; 
fileNameList = {'data_1'};
calculateTrajectoryOpt.calc2d = 0; %1 or 2 if MP/in-focus slices only
expOrSim = 'x';
clusterData = 0;
loadData = 0;
realTime = 1;

% test opt
DEBUG = []; %[1,2] for groupUnits
splitData = 0; % whether or not to split the inputData into two sets (to check the homogenity of the sample)
randomize = 0; % whether or not to randomize the order of the input data
CLUSTERMAX = 5;
downsample = 1; % 1 means no downsampling

% other
fidTxt = [];
fileList = [];

%defaults not in inputStructure (see also constants!)
writeShortFileList = 1; %writes everything on one line.
standardTags = {'spb1';'cen1'};


constants.TESTPROBABILITY = 0.95;
constants.PROB2SIDES = 1-(1-constants.TESTPROBABILITY)/2; % testprob for 2-sided test
constants.PROBOUTLIER = 0.75; % probability for ouliers. if lower, less fits get accepted
constants.PROBF = 0.7; % probability for whether linear fit is better than pause. if lower, there are less pauses
constants.STRATEGY = 1; % fitting strategy
constants.MINLENGTH = 2; % minimum length of a unit
constants.MAXDELETED = 0; % max number of deleted frames between two timepoints that is accepted
constants.MINTRAJECTORYLENGTH = 3;
constants.DEBUG = DEBUG; 
constants.DOCLUSTER = clusterData; % this does not really belong here, but it's easiest to pass it with the constants
constants.CLUSTERMIN = 1; % min # of cluster the EM algorithm looks for
constants.CLUSTERMAX = CLUSTERMAX; % max # of cluster the EM algorithm looks for
constants.CLUSTERIND = 0; % try all individual k's?
constants.CLUSTERTRY = 10; % how many times the clustering algorithm is repeated
constants.CLUSTERMINWEIGHT = 0.05; % min weight to become a significant cluster

%build list of possible identifiers
%HOME/BIODATA/SIMDATA (also: NONE/NOFILE)
%every row is: identifier, path, lengthOfPath
identifierCell = cell(3,3);

identifierCell{3,1} = 'HOME';
h = (getenv('HOME'));
if ~isempty(h) && strcmp(h(end),filesep) %should actually never happen
    h = h(1:end-1);
end
identifierCell{3,2} = h;
identifierCell{3,3} = length(identifierCell{3,2});

identifierCell{2,1} = 'BIODATA';
b = (getenv('BIODATA'));
if ~isempty(b) && strcmp(b(end),filesep)
    b = b(1:end-1);
end
identifierCell{2,2} = b;
identifierCell{2,3} = length(identifierCell{2,2});

identifierCell{1,1} = 'SIMDATA';
s = (getenv('SIMDATA'));
if ~isempty(s) && strcmp(s(end),filesep)
    s = s(1:end-1);
end
identifierCell{1,2} = s;
identifierCell{1,3} = length(identifierCell{1,2});

clear h b s
%===================================================================
%----------END DEFINE CONSTANTS ET AL--
%===================================================================




%===================================================================
%--------------TEST INPUT--------
%===================================================================


%lookfor inputData
if nargin == 0 || isempty(inputData)
    loadDataForced = 1; % we override the default
    inputData = [];
    inputDataRmIdx = [];
    inputDataIsRun = -1; % -1 : loaded input data / 0: no-run input data / 1: run input data
else
    loadDataForced = 0;
    % test inputData
    [inputData,inputDataRmIdx,inputDataIsRun] = trajectoryAnalysisIsGoodData(inputData,constants);
    
end

%go through options
%check every field for existence and change defaults - we do not need
%to consider all cases!! (whenever there is a "wrong" input, the default is taken)


if nargin < 2 || isempty(ioOpt)
    %do nothing
else
    if isfield(ioOpt,'details')
        if ioOpt.details == 1
            showDetails = 1;
        elseif ioOpt.details == 0 %assign anyway: defaults can change!
            showDetails = 0;
        end
    end
    if isfield(ioOpt,'convergence')
        if ioOpt.convergence == 1
            doConvergence = 1;
        elseif ioOpt.convergence == 0
            doConvergence = 0;
        end
    end
    if isfield(ioOpt,'saveTxt')
        if ioOpt.saveTxt == 0
            saveTxt = 0;
        elseif ioOpt.saveTxt == 1
            saveTxt = 1;
        end
    end
    
    %only now check for path
    if saveTxt && isfield(ioOpt,'saveTxtPath')
        %check whether we have been given a name or a path
        pathNameLength = length(ioOpt.saveTxtPath);
        fileSepList = strfind(ioOpt.saveTxtPath,filesep);
        if isempty(fileSepList)
            %assume we're just being given the name
            saveTxtPath = pwd;
            saveTxtName = ioOpt.saveTxtPath;
        elseif fileSepList(end)==pathNameLength
            %path ends with filesep, so no name
            saveTxtPath = ioOpt.saveTxtPath;
            saveTxtName = '';
        else
            %1:lastFileSep-1 = path, rest = filename
            saveTxtPath = ioOpt.saveTxtPath(1:fileSepList(end)-1);
            saveTxtName = ioOpt.saveTxtPath(fileSepList(end)+1:end);
        end
    end
    
    if isfield(ioOpt,'saveMat')
        if ioOpt.saveMat == 0
            saveMat = 0;
        elseif ioOpt.saveMat == 1
            saveMat = 1;
        end
    end
    
    %only now check for path
    if saveMat && isfield(ioOpt,'saveMatPath')
        %check whether we have been given a name or a path
        pathNameLength = length(ioOpt.saveMatPath);
        fileSepList = strfind(ioOpt.saveMatPath,filesep);
        if isempty(fileSepList)
            %assume we're just being given the name
            saveMatPath = pwd;
            saveMatName = ioOpt.saveMatPath;
        elseif fileSepList(end)==pathNameLength
            %path ends with filesep, so no name
            saveMatPath = ioOpt.saveMatPath;
            saveMatName = '';
        else
            %1:lastFileSep-1 = path, rest = filename
            saveMatPath = ioOpt.saveMatPath(1:fileSepList(end)-1);
            saveMatName = ioOpt.saveMatPath(fileSepList(end)+1:end);
        end
    end
    
    if isfield(ioOpt,'verbose')
        if any(ioOpt.verbose == 0) %works for [], too
            verbose = 0;
        else %if the user set something wrong, it won't have an effect
            verbose = ioOpt.verbose;
        end
    end
    
    % we do not care about the fileNameList if the input is given as run
    if isfield(ioOpt,'fileNameList') && ~inputDataIsRun
        if iscell(ioOpt.fileNameList)
            fileNameList = ioOpt.fileNameList;
            % make sure we remove fileNames if we have removed inputData
            if ~isempty(inputDataRmIdx)
                if length(fileNameList)>1
                    fileNameList(inputDataRmIdx) = [];
                else
                    %there was an error already because we killed all the
                    %inputData
                end
            end
        else
            error('ioOpt.fileNameList has to be a cell array of strings!')
        end
    end
    
    if isfield(ioOpt,'expOrSim')
        if length(ioOpt.expOrSim)>1 || ~ischar(ioOpt.expOrSim) || isempty(strfind('esx',ioOpt.expOrSim))
            error('expOrSim has to be one letter (e, s, or x)');
        else
            expOrSim=ioOpt.expOrSim;
        end
    end
    if isfield(ioOpt,'calc2d')
        calculateTrajectoryOpt.calc2d = ioOpt.calc2d;
    end
    if isfield(ioOpt,'clusterData')
        constants.DOCLUSTER = ioOpt.clusterData;
    end
    if isfield(ioOpt,'loadData')
        loadData = ioOpt.loadData;
        % if isrun is zero (not []), we force loadData to 1 if it was not 0
        % before
        if loadData && inputDataIsRun == 0
            loadData = 1;
        end
    end
    if isfield(ioOpt,'realTime')
        realTime = ioOpt.realTime;
    end
end

% == decide load data. If one of them is 1, load.
loadData = max(loadData,loadDataForced);


%check testOpt

if nargin < 3 || isempty(testOpt)
    %do nothing
else
    if isfield(testOpt,'splitData')
        splitData = testOpt.splitData;
    end
    if isfield(testOpt,'debug')
        DEBUG = testOpt.debug;
        constants.DEBUG = DEBUG;
    end
    if isfield(testOpt,'randomize')
        randomize = testOpt.randomize;
    end
    if isfield(testOpt,'clustermax')
        constants.CLUSTERMAX = testOpt.clustermax;
    end
    if isfield(testOpt,'downsample')
        downsample = max(testOpt.downsample,1); % has to be at least 1
    end
end

%===================================================================
%----------END TEST INPUT--------
%===================================================================






%===================================================================
%--------------LOAD DATA && ASK FOR PATHS---------
%===================================================================

oldDir = pwd; %remember old dir before loading

loadCt = 0;
if loadData    
    %go to biodata. remember oldDir first   
    cdBiodata(2);
    
    %buld a list of data files   
    for iInput = 1:loadData
        fileList = loadFileList({'*.mte;*.mts;*.mtx','results files';...
                '*-data-??-???-????-??-??-??.mat','project data files'},[2;1],...
                [],['load data for run - ' num2str(loadCt+1) ' -']);
            if ~isempty(fileList)
                % only count loaded data if we have loaded at all
                loadCt = loadCt + 1;
                inputLoaded(loadCt).fileList = fileList;
            end    
    end
else
    inputLoaded = [];
end

% check whether we do have data at all
if loadCt == 0 && isempty(inputData)
    disp('no files loaded - end evaluation')
    return
end

% count runs (we do not count the additional runs yet!)
if inputDataIsRun == 0
    numRun = 1;
else
    numRun = loadCt + length(inputData);
end


%---interrupt load data and ask for paths - calculating the trajectories
%takes too much time!

%first: textFile - only if we have just one run!
if saveTxt && numRun == 1
    if isempty(saveTxtName) %only ask if no name has been given before
        helpTxt = 'Please select filename,';
        if isempty(saveTxtPath)
            helpTxt = [helpTxt, ' pathname,'];
        else
            cd(saveTxtPath)
        end
        helpTxt = [helpTxt, ' and filetype to save the summary of the results as text file! If you press ''cancel'', the file will not be saved'];
        
        %tell the user what's going on
        ans = myQuestdlg(helpTxt,'','OK','cancel','OK');
        if strcmp(ans,'OK')
            cdBiodata(0);
            [saveTxtName,saveTxtPath,eos] = uiputfile({'*.mte','experimental MT data';...
                    '*.mts','simulation MT data';...
                    '*.mtx','any (mixed) MT data'},'save results as text file');
        else
            saveTxtName = 0;
        end
        
        %if user cancelled, nothing will be save                               
        if saveTxtName == 0
            saveTxt = 0;
        else
            
            %get expOrSim
            eosVect = 'esx';
            expOrSim = eosVect(eos);
        end
    end
    
    %append the file extension
    if saveTxt && ~strcmp(saveTxtName(end-3),'.')
        saveTxtName = [saveTxtName,'.mt',expOrSim];
    else
        %the user has chose something him/herself
    end
else
    % we do not save at all
    saveTxt = 0;
    saveTxtName = 0;
end

%second: matFile - save anyway
if saveMat
    if isempty(saveMatName) %only ask if no name has been given before
        helpTxt = 'Please select filename';
        if isempty(saveMatPath)
            helpTxt = [helpTxt, ' and pathname,'];
        else
            cd(saveMatPath)
        end
        helpTxt = [helpTxt, ' to save the results as mat-file! If you press ''cancel'', the file will not be saved'];
        
        %tell the user what's going on
        ans = myQuestdlg(helpTxt,'','OK','cancel','OK');
        if strcmp(ans,'OK')
            
            [saveMatName,saveMatPath] = uiputfile({'*.mat;*.mt*','MT-dynamics files'},'save results as mat file');
        else
            saveMatName = 0;
        end
        
        %if user cancelled, nothing will be save                               
        if saveMatName == 0
            saveMat = 0;
        end
        
    end
end

%resume loading data - downsampling is done in tALoadDAta
   


%========= LOAD DATA && SET UP RUNS ===============================
% new data is appended to old one
loadOptions = struct('standardTags', {standardTags},...
    'downsample', downsample,...
    'splitData', splitData,...
    'randomize', randomize,...
    'calculateTrajectoryOpt',calculateTrajectoryOpt,...
    'realTime',realTime);

[run, fileNameListSave] = trajectoryAnalysisLoadData(...
    inputLoaded, constants, inputData, fileNameList, loadOptions);
%===================================================

cd(oldDir);

clear fileList helpTxt problem i inputData fileNameList loadOptions

%===================================================================
%----------END LOAD DATA  && ASK FOR PATHS---------
%===================================================================



%===================================================================
%--------------WRITE FILE-LIST--------
%===================================================================

%already write the fileList for the results-file (if selected): if
%something happens during the calculations, at least the list is not lost.
if saveTxt
    %create the file
    fidTxt = fopen([saveTxtPath,saveTxtName],'w');
    
    %if we selected to write the short version of the fileList, we do not
    %use a line break between identifier/options and the rest of the
    %filename
    if writeShortFileList
        separationString = '   ';
    else
        separationString = '\n';
    end
    
    %write introduction to file
    fprintf(fidTxt,'%s\n%s\n%s\n%s\n','%  MICROTUBULE DYNAMICS ANALYSIS - list of filenames',...
        '%    this file contains a list of filenames that can be used for MT dynamic analysis',...
        '%    the list is made up as: Identifier#{''tag1'',''tag2''}# \n  restOfPathIncludingFileName, where Identifier is the environment variable for the particular path',...
        '%    Please do not uncomment this header or delete the ''***'' that mark the end of the filenames. Fileseps can be windows or linux type.');
    
    
    %now loop through the fileNameList, find the identifier and write the
    %file
    for nFile = 1:length(fileNameListSave)
        
        %init variables
        identifier = '';
        restOfFileName = '';
        tagList = {''};
        
        %read tagList
        tagList = run(1).data(nFile).info.tags;
        
        %read fileName
        fileName = fileNameListSave{nFile};
        lengthFileName = length(fileName);
        
        %now sieve the fileName until we find the identifier, or know for
        %sure that there isn't any
        %use if/elseif, because we want biodata/simdata to override home
        if ~isempty(identifierCell{1,2}) && strcmpi(identifierCell{1,2},fileName(1:min(identifierCell{1,3},lengthFileName))) %check for SIMDATA
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{1,1};
            restOfFileName = fileName(identifierCell{1,3}+1:end);
            
        elseif ~isempty(identifierCell{2,2}) && strcmpi(identifierCell{2,2},fileName(1:min(identifierCell{2,3},lengthFileName))) %check for BIODATA
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{2,1};
            restOfFileName = fileName(identifierCell{2,3}+1:end);
            
        elseif ~isempty(identifierCell{3,2}) && strcmpi(identifierCell{3,2},fileName(1:min(identifierCell{3,3},lengthFileName))) %check for HOME
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{3,1};
            restOfFileName = fileName(identifierCell{3,3}+1:end);
            
        elseif exist(fileName,'file') %check for NONE
            
            %assign none
            identifier = 'NONE';
            restOfFileName = fileName;
            
        else %we have no valid filename at all
            
            identifier = 'NOFILE';
            restOfFileName = fileName;
            
        end %check for identifier and restOfFilename
        
        %now write everything to file
        fprintf(fidTxt,['%s#%s#%s',separationString,'%s\n'],...
            identifier, tagList{1}, tagList{2}, restOfFileName);
        
        
    end %for nFile = 1:length(fileNameList)
    
    %close off the file by writing '***'
    fprintf(fidTxt,'\n***\n\n');
    
    %do not close file yet - we will write more!
end
%===================================================================
%----------END WRITE FILE-LIST--------
%===================================================================




%===================================================================
%--------------CALCULATE---------
%===================================================================

if ~isempty(DEBUG) || ~saveTxt % if no file, we do not care
    for iRun = 1:length(run)
        trajectoryDescription(iRun,1) =...
            trajectoryAnalysisMain(run(iRun).data,...
            constants,showDetails,doConvergence,verbose,run(iRun).fileNameList);
    end
else
    try
        for iRun = 1:length(run)
            trajectoryDescription(iRun,1) =...
                trajectoryAnalysisMain(run(iRun).data,...
                constants,showDetails,doConvergence,verbose,run(iRun).fileNameList);
        end
    catch
        if ~isempty(fidTxt) 
            fclose(fidTxt);
        end
        rethrow(lasterr)
    end
end
%===================================================================
%----------END CALCULATE---------
%===================================================================


%-----add additional info for trajDes here
%
%----------------------------------------

%===================================================================
%-------------STORE DATA---------
%===================================================================

%splitData only works for saveMat

%save mat file first, because it's less lines
if saveMat
    save([saveMatPath,saveMatName],'trajectoryDescription');
end

%save text file
if saveTxt
    
    for iRun = 1:length(run)
        
        %get fieldNames for saving
        statisticsTitles = fieldnames(trajectoryDescription(iRun).individualStatistics(1).summary);
        numStats = length(statisticsTitles);
        
        %---write overall statistics, date, probabilities
        fprintf(fidTxt,'\n\n---OVERALL STATISTICS, run %i or %i---\n', iRun, length(run));
        fprintf(fidTxt,[nowString,'\n']);
        fprintf(fidTxt,'Slope   Test: %1.3f\n',constants.TESTPROBABILITY);
        fprintf(fidTxt,'Pause   Test: %1.3f\n',constants.PROBF);
        fprintf(fidTxt,'Outlier Test: %1.3f\n\n\n',constants.PROBOUTLIER);
        %fprintf(fidTxt,'Clustering:')
        
        %read them first
        if doConvergence
            overallStats = trajectoryDescription(iRun).convergenceStatistics(end);
        else
            overallStats = trajectoryDescription(iRun).overallStatistics;
        end
        %make cell for writing down
        statisticsCell = [statisticsTitles,struct2cell(overallStats)];
        
        for nStat = 1:numStats
            %(I so love these formatted strings)
            txt2write  = statisticsCell{nStat,1};
            vars2write = [statisticsCell{nStat,2:end}];
            formatString = '%25sVARIABLES\n';
            % allow for any nunber of vars2write
            formatString = strrep(formatString,'VARIABLES',repmat('\t%5.6f',[1,length(vars2write)]));
            fprintf(fidTxt,formatString,txt2write,vars2write);
        end
        
        %---loop through the individual trajectories, print their overall statistics
        fprintf(fidTxt,'\n\n---individual data---');
        numFiles = length(run(iRun).fileNameList);
        for nFile = 1:numFiles
            %write a nice title
            fprintf(fidTxt,'\n\noverall statistics for %s\n',run(iRun).fileNameList{nFile});
            
            %read the list
            statisticsCell = [statisticsTitles,struct2cell(trajectoryDescription(iRun).individualStatistics(nFile).summary)];
            for nStat = 1:numStats
                txt2write  = statisticsCell{nStat,1};
                vars2write = [statisticsCell{nStat,2:end}];
                formatString = '%25sVARIABLES\n';
                formatString = strrep(formatString,'VARIABLES',repmat('\t%5.6f',[1,length(vars2write)]));
                fprintf(fidTxt,formatString,txt2write,vars2write);
            end
        end
        %DO NOT UNCOMMENT - USE ONLY AS TEMPLATE           
        %         %loop through individual trajectories, print measurements
        %         fprintf(fidData,'\n\n---individual data II---');
        %         for i = 1:numData
        %             fprintf(fidData,'\n\ndetailed statistics for %s\n',mtdDat(i).info.name);
        %             
        %             
        %             %group
        % %             numEntries = size(mtdDat(i).stateList.group,1);
        % %             fprintf(fidData,'group data [counter, state(1/2/-1/-2/0), startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2, avg. time in undetermined state, sum of single counters, significance1, significance2]\n');
        % %             for j = 1:numEntries
        % %                 fprintf(fidData,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n*',mtdDat(i).stateList.group(j,:));
        % %             end
        %         end
        %-----------------------
    end % for iRun = 1:length(run)
    fclose(fidTxt);
end %if saveTxt


%===================================================================
%---------END STORE DATA---------
%===================================================================
