function run = loadRunsFromFile(nRunsOrFileList,addOri,addPos,addSig0,addLag,addDist,addTime,addTP)
%LOADRUNSFROMFILE lets the user load runs from file
%
% SYNOPSIS run = loadRunsFromFile(nRunsOrFileList,addOri,addPos,addSigma0,addTime)
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
%          add...          (all opt) Switches to add more data to the
%                            output: orientation, position 
%                            structure, sigma0, timeLapse (all default 0),
%                            distance, time, timePoints (all default 1)
%                            
%
% OUTPUT   run             input structure for trajectoryAnalysis with e.v.
%                            additional fields
%
% c: jonas, 05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================
% TEST INPUT
%================

% SWITCHES
if nargin < 8 | isempty(addTP)
    addTP = 1;
end
if nargin < 7 | isempty(addTime)
    addTime = 1;
end
if nargin < 6 | isempty(addDist)
    addDist = 1;
end
if nargin < 5 | isempty(addLag)
    addLag = 0;
end
if nargin < 4 | isempty(addSig0)
    addSig0 = 0;
end
if nargin < 3 | isempty(addPos)
    addPos = 0;
end
if nargin < 2 | isempty(addOri)
    addOri = 0;
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
            [tmpData,ori,pos,sig0,tL] = calculateTrajectoryFromIdlist(idlist2use,allDat.dataProperties,tag1,tag2,calculateTrajectoryOpt);
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
            if addLag
                data(dataCt).timeLapse = tL;
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
        runCt = runCt +1;
        dataCt = 1;
    end
    
end %  for iRun = runCt:length(fileListStruct)

    

    