function run = loadRunsFromFile(nRunsOrFileList,varargin)
%LOADRUNSFROMFILE lets the user load runs from file
%
% SYNOPSIS run = loadRunsFromFile(nRunsOrFileList,varargin)
%
% INPUT    nRunsOrFilelist (opt) if scalar, then n Runs are loaded.
%                            otherwise, specify inputList a 1-by-nRuns
%                            structure with field fileList, which is a
%                            struct with fields .file, .type of length n, 
%                            where n is the number of files to be loaded. 
%                               .file  : full path for the file including filename
%                               .type  : (numeric) filetype (results files: 1, projFiles: 2).
%                               .opt   : cell array with additional information, as
%                                        stored on the identifier-line. opt{1} is the
%                                        identifier
%                            Default : 1
%
%          varargin        (all opt) Strings to select which data to add to the output.
%                             If the string starts with a '-', the data is
%                             not loaded. If it starts with nothing, or
%                             '+', the data is added.
%
%                             Strings   Explanation     Default
%                               dist     distance           1
%                               time     time [s]           1
%                               tp       timePoints         1
%                               ori      orientation        0
%                               sig0     sigmaZero          0
%                               dp       dataProperties     0
%                               SNR      SNR max            0
%                               isT      tracked or not     0
%                               pos      position structure 0
%                            
%
% OUTPUT   run             input structure for trajectoryAnalysis with e.v.
%                            additional fields. Fieldnames:
%                      
%                       .fileNameList
%                       .data
%                           .distance
%                           .time
%                           .timePoints
%                           .orientation
%                           .sigmaZero
%                           .dataProperties
%                           .snrMax
%                           .isTracked
%                           .position
%
% c: jonas, 05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================
% TEST INPUT
%================

% defaults
addTP   = 1;
addTime = 1;
addDist = 1;
addDP   = 0;
addSig0 = 0;
addPos  = 0;
addOri  = 0;
addSnr  = 0;
addIsT  = 0;

if nargin == 0 || isempty(nRunsOrFileList)
    error('Specify nRunsOrFileList in loadRunsFromFile!')
end

if nRunsOrFileList == 0
    return
end


% SWITCHES
for in = 1:length(varargin)
    arg2check = varargin{in};
    
    % make sure it's a string
    % FUTURE: if it's the last argument, check whether it is some option
    % structure
    if ~isstr(arg2check)
        error('Options for loadRunsFromFile have to be strings! (offending argument#: %i)',in);
    end
    
    % check for + or - sign
    switch arg2check(1)
        case '+'
            newValue = 1;
            arg2check = arg2check(2:end);
        case '-'
            newValue = 0;
            arg2check = arg2check(2:end);
        otherwise
            newValue = 1;
    end
    
    % find which to set
    switch arg2check
        case 'dist'
            % distance 
            addDist = newValue;
        case 'time'     
            % time [s] 
            addTime = newValue;
        case 'tp'
            % timePoints 
            addTP = newValue;
        case 'ori'
            % orientation    
            addOri = newValue;
        case 'sig0'
            % sigmaZero
            addSig0 = newValue;
        case 'dp'
            % dataProperties
            addDP = newValue;
        case 'SNR'
            % SNR max
            addSnr = newValue;
        case 'isT'
            % tracked or not
            addIsT = newValue;
        case 'pos'
            % add position structure
            addPos = newValue;
         
        otherwise
            warning('Option %i for loadRunsFromFile not recognized',in);
    end
end

% % tags
% if nargin < 2 | isempty(tags)
%     tags = {'spb1';'cen1'};
% elseif iscell(tags) & length(tags)==2
%     % good
% else
%     error('Please specify two tags in a cell')
% end
standardTags = {'spb1';'cen1'};


%nRunsOrFilelist
if nargin==0 | isempty(nRunsOrFileList)
    loadData = 1;
    nRuns    = 1;
elseif isnumeric(nRunsOrFileList) & isfinite(nRunsOrFileList)
    loadData = nRunsOrFileList;
    nRuns    = nRunsOrFileList;
elseif isstruct(nRunsOrFileList) & isfield(nRunsOrFileList,'fileList')
    loadData = 0;
    inputList = nRunsOrFileList;
    nRuns = length(nRunsOrFileList);
else
    error('unable to handle first input argument')
end

clear nRunsOrFileList;

%================

%================
% load fileLists
%================

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
            inputList(loadCt).fileList = fileList;
        end    
    end
    
    % check whether we do have data at all
    if loadCt == 0 
        disp('no files loaded - end evaluation')
        return
    end
    
    % update nRuns
    nRuns = loadCt;
    
end

cd(oldDir);
%================



%================
% LOAD FILES
%================
problem = [];
dataCt = 1;
runCt = 1;
for iRun = 1:nRuns
    fileList = inputList(iRun).fileList;
    for iFile = 1:length(fileList)
        
        
        try
            %load all
            allDat = load(fileList(iFile).file);
            
            %load the idlist specified in lastResult
            eval(['idlist2use = allDat.',allDat.lastResult,';']);
            
            % if addSnr, we add the idlist_L, if it exists
            if addSnr
                if isfield(allDat,'idlist_L')
                    calculateTrajectoryOpt.oldIdlist = allDat.idlist_L;
                end
            end
                    
            
            
            %---prepare calculate trajectory
            
            %choose tags
            %check whether there is a non-empty field opt, else just use
            %standardTags
            
            flopt = fileList(iFile).opt; %increase readability
            if isempty(flopt) | length(flopt)<3  | isempty(flopt{2}) | isempty(flopt{3})%there could be more options in the future!
                % use default tag1, tag2
                tag1 = standardTags{1};
                tag2 = standardTags{2};
            else
                tag1 = fileList(iFile).opt{2};
                tag2 = fileList(iFile).opt{3};
            end
            
            %store identifier
            if isempty(fileList(iFile).opt) | isempty(fileList(iFile).opt{1})
                calculateTrajectoryOpt.identifier = 'NONE';
            else
                calculateTrajectoryOpt.identifier = fileList(iFile).opt{1};
            end
            
            %-----calculate trajectory -- the assignment data(i) = output.a/b/c does not work if data is []!!
            [tmpData,ori,pos,sig0,dataProperties,snrMax,isTracked] = calculateTrajectoryFromIdlist(idlist2use,allDat.dataProperties,tag1,tag2,calculateTrajectoryOpt);
            %-------------------------
            
            % add standard data
            if addDist
                data(dataCt).distance = tmpData.distance;
            end
            if addTime
                data(dataCt).time     = tmpData.time;
            end
            if addTP
                data(dataCt).timePoints = tmpData.timePoints;
            end
            
            % add more data
            if addOri
                data(dataCt).orientation = ori;
            end
            if addPos
                data(dataCt).position = pos;
            end
            if addSig0
                data(dataCt).sigmaZero = sig0;
            end
            if addDP
                data(dataCt).dataProperties = dataProperties;
            end
            if addSnr
                data(dataCt).snrMax = snrMax;
            end
            if addIsT
                data(dataCt).isTracked = isTracked;
            end
            
            %remember fileName
            fileNameList{dataCt} = fileList(iFile).file;
            
            % prepare runCt for next turn
            dataCt = dataCt + 1;
            
            
        catch
            if isempty(problem)
                problem = lasterr;
            end
            disp([fileList(iFile).file, ' could not be loaded:',char(10),problem])
        end
        
        clear lr id dp
    end %for i = 1:length(fileList)
    
    % make sure we loaded something
    if dataCt == 1
        error('no data loaded!') % alternatively, we just could move on without updating runCt/dataCt
    else
        run(runCt).data = data;
        run(runCt).fileNameList = fileNameList;
        clear data
        clear fileNameList;
        runCt = runCt +1;
        dataCt = 1;
    end
    
end %  for iRun = runCt:length(fileListStruct)

    

    