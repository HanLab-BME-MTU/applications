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
%                               nanDist  distance as nanList
%                                 output will have field
%                                 timeInterval instead of
%                                 time/timePoints. Default: 0
%                               idlist   entire idlist      0
%                               realTime use real timestamp 1
%
%                              In order to supply a condition idlists have
%                              to fulfill, please input a structure array
%                              as any of the variable input arguments.
%                              Fields:
%                                .check
%                                .askOptions
%                              see checkIdlist for the options for these
%                              arguments.
%
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
%                           .idlist
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
convDist= 0;
addIdlist = 0;
useRealTime = 1;
checkStructIn = struct('check',[],'askOptions',[]);


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
    if ~ischar(arg2check)
        if isstruct(arg2check)
            if isfield(arg2check,'check')
                checkStructIn.check = arg2check.check;
            end
            if isfield(arg2check,'askOptions')
                checkStructIn.askOptions = arg2check.askOptions;
            end
        else
            error('Options for loadRunsFromFile have to be strings! (offending argument#: %i)',in);
        end
    else

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
            case 'nanDist'
                convDist = newValue;
            case 'idlist'
                addIdlist = newValue;
            case 'realTime'
                useRealTime = newValue;

            otherwise
                warning('Option %i (''%s'') for loadRunsFromFile not recognized',in,arg2check);
        end
    end
end

% check whether we have to call calculateTrajectoryFromIdlist
selectionVector = [addTP;addTime;addDist;addDP;addSig0;addPos;addOri;...
    addSnr;addIsT;convDist;addIdlist];
if all(selectionVector == [zeros(length(selectionVector)-1,1);1])
    calcTraj = 0;
else
    calcTraj = 1;
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
if nargin==0 || isempty(nRunsOrFileList)
    loadData = 1;
    nRuns    = 1;
elseif isnumeric(nRunsOrFileList) && isfinite(nRunsOrFileList)
    loadData = nRunsOrFileList;
    nRuns    = nRunsOrFileList;
elseif isstruct(nRunsOrFileList) && isfield(nRunsOrFileList,'fileList')
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
dataCt = 0;
runCt = 0;
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
            if isempty(flopt) || length(flopt)<3  || isempty(flopt{2}) || isempty(flopt{3})%there could be more options in the future!
                % use default tag1, tag2
                tag1 = standardTags{1};
                tag2 = standardTags{2};
            else
                tag1 = fileList(iFile).opt{2};
                tag2 = fileList(iFile).opt{3};
            end

            %store identifier
            if isempty(fileList(iFile).opt) || isempty(fileList(iFile).opt{1})
                calculateTrajectoryOpt.identifier = 'NONE';
            else
                calculateTrajectoryOpt.identifier = fileList(iFile).opt{1};
            end

            % check whether to convert data to nanList
            if convDist
                calculateTrajectoryOpt.nanList = 1;
            end

            % realTime switch
            calculateTrajectoryOpt.realTime = useRealTime;
            
            % check idlist. Properly prepare input. We want, in any case,
            % check for tag1&tag2. Will throw an error if not n-by-3 cell
            % input for checkCell or check.
            checkStruct = struct('check',[],'askOptions',[]);
            if isempty(checkStructIn.check)
                checkStruct.check = {7,tag1,tag2};
            else
                if iscell(checkStructIn.check)
                    checkStruct.check = [checkStructIn.check;{7,tag1,tag2}];
                elseif ischar(checkStructIn.check)
                    if strcmp(checkStructIn.check,'ask')
                        checkStruct.check = 'ask';
                        if ~isfield(checkStruct.askOptions,'checkCell')
                            checkStruct.askOptions.checkCell = {7,tag1,tag2};
                        else
                            checkStruct.askOptions.checkCell = [checkStructIn.askOptions.checkCell;{7,tag1,tag2}];
                        end
                    else
                        checkStruct.check = [{checkStructIn.check,'',''};{7,tag1,tag2}];
                    end
%                 else
%                     checkStruct.check = checkStructIn.check;
                end
            end
            [goodIdlist,errorMessage,goodTimes] = checkIdlist(idlist2use,checkStruct.check,checkStruct.askOptions);
            % quit with error if necessary
            error(errorMessage);
            % remove bad goodTimes
            [idlist2use(~goodTimes).linklist] = deal([]); %#ok<AGROW>

            if calcTraj
                %-----calculate trajectory -- the assignment data(i) = output.a/b/c does not work if data is []!!
                [tmpData,ori,pos,sig0,dataProperties,snrMax,isTracked] =...
                    calculateTrajectoryFromIdlist(...
                    idlist2use,allDat.dataProperties,tag1,tag2,calculateTrajectoryOpt);
                %-------------------------
            end

            % update dataCt
            dataCt = dataCt + 1;

            % add standard data
            if addDist
                data(dataCt).distance = tmpData.distance;
            end
            if convDist
                %                 data(dataCt).timeInterval = tmpData.timeInterval;
                % something is wrong here. This is just a patch until Jonas takes a look at
                % this
                data(dataCt).timeInterval = tmpData.time;
            else
                if addTime
                    data(dataCt).time     = tmpData.time;
                end
                if addTP
                    data(dataCt).timePoints = tmpData.timePoints;
                end
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
            if addIdlist
                data(dataCt).idlist = idlist2use;
            end
            
            %KJ: add tag names in order to keep track
            data(dataCt).tags{1} = tag1;
            data(dataCt).tags{2} = tag2;

            %remember fileName
            fileNameList{dataCt} = fileList(iFile).file;




        catch
            if isempty(problem)
                problem = lasterr;
            end
            disp([fileList(iFile).file, ' could not be loaded:',char(10),problem])
        end

        clear lr id dp
    end %for i = 1:length(fileList)

    % make sure we loaded something
    if dataCt == 0
        error('no data loaded!') % alternatively, we just could move on without updating runCt/dataCt
    else
        runCt = runCt +1;
        run(runCt).data = data;
        run(runCt).fileNameList = fileNameList;
        clear data
        clear fileNameList;

        dataCt = 0;
    end
    
    % remove checkIdlist from memory to clear persistent variables
    clear checkIdlist

end %  for iRun = runCt:length(fileListStruct)


