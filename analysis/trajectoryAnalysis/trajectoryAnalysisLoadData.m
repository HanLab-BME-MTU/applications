function [run,fileNameListSave] = trajectoryAnalysisLoadData(fileListStruct,constants,optData,fileNameList,loadOptions)
%TRAJECTORYANALYSISLOADDATA loads trajectories given in fileList into data
%
% SYNOPSIS [run,fileNameListSave] = trajectoryAnalysisLoadData(fileListStruct,optData,fileNameList,loadOptions)
%
% INPUT    fileListStruct(1:nfl).fileList     : as returned from loadFileList (can be empty)
%          constants    : needs field constants.MINTRAJECTORYLENGTH
%          optData      : (opt) data already structured like the output
%                         structure (run) or single set of data
%          fileNameList : list of fileNames (cell array of strings!)
%                         only if optData is not a run structure
%          loadOptions  : (opt) data structure for additional options
%                   the first four are applied in this order
%               - splitData    [{0}/1] whether to split the data randomly
%                              into two exclusive sets (exclusive to subset)
%                              overrides subset. 
%               - subset       [fraction,n] {[1,1]} randomly takes n unique
%                              subsets of size fraction*nData. MaxN is 100
%                              (exclusive to splitData). 
%               - downsample   scalar. downsample n means that out of one
%                              trajectory, n trajectories will be made by
%                              taking every nth data point from the
%                              original trajectory 
%               - randomize    [{0},1] whether to randomize the order of the data
%
%               - standardTags {standardTag1; standardTag2} tag names given
%                              to the two points at either end of the
%                              distance
%               - calculateTrajectoryOpt : options structure for
%                                          calculateTrajectory, with fields
%                           .info :   struct that should be returned in the
%                                     output field info. 
%                           .calc2d : [{0}/1/2] depending on whether the
%                                     normal 3D-data or the maxProjection
%                                     or the in-focus-slice data should be
%                                     used.
%
%
% OUTPUT   run(1:nRun): structure with data for every individual run with
%           fields
%               data(1:n)  : structure containing n different trajectories with fields
%                   - distance   tx2 array [distance, sigmaDistance] in microns
%                   - time       tx2 array [time, sigmaTime] in seconds
%                   - timePoints tx1 array [timePoint#] 
%                   - info       struct containing additional information about the data
%                       -tags    : cell containing the two strings designating
%                                  the tags between which the distance is
%                                  measured, e.g. {'spb1','cen1'}
%               fileNameList : old filenames with suffix _ds#downsample_#nth
%
%          fileNameListSAve    : original list of files for saving data
%
% pulled out of trajectoryAnalysis 3/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set defaults
downsample = 1;
splitData  = 0;
randomize  = 0;
subset     = [1 1];
standardTag1 = 'spb1';
standardTag2 = 'cen1';
standardFileName = {'data_1'};


%============
% test input
%============

if nargin < 1 | (isempty(fileListStruct) & (nargin < 3 | isempty(optData)))
    error('not enough input arguments to load data')
end

% we do not test fileList here - if it is crap, nothing will be loaded (hopefully)

if nargin < 2 | isempty(constants)
    error('please specify constants.MINTRAJECTORYLENGTH')
end



% test loadOptions before testing optData - there could be standardTags
if nargin < 5 | isempty(loadOptions)
    % set all to standard
else
    if isfield(loadOptions,'downsample')
        % downsample is at least 1
        downsample = max(loadOptions.downsample,1);
    end
    
    if isfield(loadOptions,'splitData')
        splitData = loadOptions.splitData;
    end
    
    if isfield(loadOptions,'subset')
        subset = loadOptions.subset;
        if subset(1) <= 0 | subset(1) > 1
            error('please input the first parameter of loadOptions.subset as a fraction between 0+ and 1')
        end
        if length(subset) == 1
            error('please specify how many subsets you want to have');
        end
        if subset(2)> 100
            error('You can''t choose more than 100 subsets at a time, sorry.')
        end
    end
    if isfield(loadOptions,'standardTags')
        if length(loadOptions.standardTags) == 2 & iscell(loadOptions.standardTags)
            standardTag1 = loadOptions.standardTags{1};
            standardTag2 = loadOptions.standardTags{2};
        end
    end
    if isfield(loadOptions,'calculateTrajectoryOpt')
        calculateTrajectoryOpt = loadOptions.calculateTrajectoryOpt;
    end
    if isfield(loadOptions,'randomize')
        randomize = loadOptions.randomize;
    end
        
end

rmOptDataIdx = [];
if nargin < 3 | isempty(optData)
    % init optData
    optData = [];
    optDataIsRun = -1;
else
    % test optData
    try
        [optData,rmOptDataIdx,optDataIsRun] = trajectoryAnalysisIsGoodData(optData,constants,standardTag1,standardTag2);
    catch
        % if there was just not enough data, but we are loading, it's not
        % so much of a problem
        if strcmp(lasterr,'no valid data points after removal')
            if ~isempty(fileList)
                warning('TRAJECTORYANALYSIS:inputTest','no data left in optData after removal')
                optData = [];
                fileNameList = []; % init or kill
            else
                % if there is no data, we are not happy
                rethrow(lasterror)
            end
        else
            rethrow(lasterror)
        end
    end
end

% test fileNameList after testing optData - it could be a run!
if nargin < 4
    fileNameList = []; % just init
else
    if  ~isempty(optData) & optDataIsRun == 0
        
        %set default if necessary
        if isempty(fileNameList)
            fileNameList = standardFileName;
        else
            if ~iscell(fileNameList) | ~isstr(fileNameList{1})
                error('fileNameList has to be a cell array of strings!')
            end
        end
        
        % make sure we have the right number of filenames
        
        lengthFNL = length(fileNameList);
        lengthData = length(optData);
        if (lengthData > lengthFNL) 
            %if more data than fnl : fill fnl if length 1
            %if > 1 return error
            %if less data than fnl : return warning
            if lengthFNL == 1
                if isempty(fileNameList{1}) || strcmp(fileNameList{1},'data_1');
                    dataStr = [repmat('data_',[lengthData,1]),num2str([1:lengthData]')]; %string data_  1 - data_999
                    dataStr = regexprep(dataStr,' ','_'); %dataStr = data___1 - data_999
                    fileNameList = cellstr(dataStr);
                else
                    fileNameList = repmat(fileNameList,[lengthData,1]);
                end
            else
                error('not enough fileNames in ioOpt.fileNameList!')
            end
        elseif lengthData < lengthFNL
            % maybe we have just removed some of the optData
            if lengthFNL == lengthData + length(rmOptDataIdx)
                fileNameList(rmOptDataIdx) = [];
            else
                warning('there are more fileNames than data - assignment of names could be inaccurate!')
            end
        end
    end
end


%=========================================================================


%=======================================================
% LOAD DATA
%=======================================================

% append to optData if it is a trajectoryData, otherwise, append to the
% runs

switch optDataIsRun
    case 0 % trajectoryData is passed down
        
        if isempty(fileListStruct)
            % we do not load, so just write run
            run = struct('data',optData,'fileNameList',{fileNameList});
            runCt = 1;
            dataCt = 1;
        else
            
            run = struct('data',[],'fileNameList',[]);
            data = optData; % init data only here, we can not init it empty
            runCt = 1;
            dataCt = length(data) + 1; % we append data
        end
        
    case 1 % run is passed down
        
        run = optData;
        runCt = length(run) + 1; % we append runs
        dataCt = 1;
        
    case -1 % nothing is passed down
        
        run = struct('data',[],'fileNameList',[]);
        runCt = 1;
        dataCt = 1;
        
end


problem = [];
%load the data
for iRun = 1:length(fileListStruct)
    fileList = fileListStruct(iRun).fileList;
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
                tag1 = standardTag1;
                tag2 = standardTag2;
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
            data(dataCt) = calculateTrajectoryFromIdlist(idlist2use,allDat.dataProperties,tag1,tag2,calculateTrajectoryOpt);
            %-------------------------
            
            % test whether we want to keep this
            if length(data(dataCt).timePoints)<constants.MINTRAJECTORYLENGTH
                error('not enough data points')
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

% remember the filenames for all the data loaded before we tinker with the
% data
if length(run) == 1
    fileNameListSave = run(1).fileNameList;
else
    fileNameListSave = {}; % we do not need this if we do not save
end

%===================================================================
%----------SPLIT DATA SETS / SUBSET
%===================================================================
if splitData
    nRuns = length(run);
    % make indexLists to store runs
    oddRuns = [1:nRun]*2-1;
    evenRuns = oddRuns +1;
    for iRun = 1:nRuns
        nData = length(run(iRun).data);
        if mod(nData,2)~=0
            disp('warning: odd number of files, not taking into account last file')
            nData = nData - 1;
        end
        % create two sets (thanks to Aaron for the idea!)
        % get nData random numbers, sort, take the first nData/2 indices for
        % the first half etc.
        rankList = randperm(nData);
        dataSet1 = rankList(1:nData/2);
        dataSet2 = rankList(nData/2+1:end);
        
        newRun(oddRuns(iRun)).data  = run(iRun).data(dataSet2);
        newRun(evenRuns(iRun)).data  = run(iRun).data(dataSet1); 
        newRun(oddRuns(iRun)).fileNameList = run(iRun).fileNameList(dataSet2);
        newRun(evenRuns(iRun)).fileNameList  = run(iRun).fileNameList(dataSet1); % same as above
        
    end % for iRun = 1:nRuns
elseif ~all(subset == [1 1])
    % subset = [frac,n]. We take n subsets of size frac*nData, or the
    % maximum number of possible subsets for the given fraction. 
    % max # of subsets: 
    % size subset s=round(frac*nData) or 1 or nData
    % max# = nData!/(s! * (nData-s)!)
    %
    % we loop through the runs, and fill newRuns. We init first with n
    % subsets per run, but might have to throw away some entries at the end
    %
    % to select the random subsets, we first make a lookup-table. Then we
    % randomly select one of the table entries and delete this entry.
    % Rinse, repeat. Of course, it's much faster to throw away entries in
    % the list if n is close to the maximum. 
    % OK, this is a bad idea. The number of possible permutations is just
    % too high to create the LUT. Instead, we will limit the number of
    % allowed permutations to 100, make 200 random permutations, and select
    % the first n unique ones. This is quite quick and works probably
    % always, except if you have a very bad day.
    
    % turn off warnings
    oldWarnings = warnings;
    warning off MATLAB:nchoosek:LargeCoefficient
    
    % preassign newRun
    nRuns = length(run);
    newRun(1:nRuns*subset(1)) = struct('data',[],'fileNameList',[]);
    
    newRunCt = 0;
    % loop through runs
    for iRun = 1:nRuns
        
        % calculate how many how large subsets we want to make. Minimum size: 1
        nData = length(run(iRun).data);
        subSetSize = max(round(subset(1)*nData),1);
        numSubSets = subset(2);
        
        % calculate how many subsets are possible
        maxNumSS = nchoosek(nData,subSetSize);
        
        %...and adjust numSubSets if necessary
        numSubSets = min(maxNumSS,numSubSets);
        
        enoughSubsets = 0;
        numAttempts = 0;
        maxNumAttempts = 3;
        
        while enoughSubsets & numAttempts<maxNumAttempts
            
            numAttempts = numAttempts + 1;
            
            % now create 200 randperms
            randomSets = zeros(200,nData);
            for i = 1:200
                randomSets(i,:)=randperm(nData);
            end
            
            % select the subset
            poSS = randomSets(:,1:subSetSize);
            poSS = sort(poSS,2);
            poSS = unique(poSS,'rows');
            % now we have the possible subset. Check if we have enough subsets.
            numPossibleSS = size(poSS,1);
            enoughSubsets = >=numSubSets;
            
        end
        
        % we tried now up to three times, so if we still do not have enough
        % subsets, that's just pure bad luck
        
        if ~enoughSubsets
            warning('sorry, not enough subsets could be found. Using all we have')
            selectedSS = poSS;
            numSubSets = numPossibleSS;
        else
            % careful! unique sorts
            rp = randperm(numPossibleSS);
            selectedSS = poSS(rp(1:numSubSets),:);
        end
        
        % loop through all the subsets and assign
        for iSet = 1:numSubSets
            newRunCt = newRunCt + 1;
            newRun(newRunCt).data = run(iRun).data(selectedSS(iSet,:));
            newRun(newRunCt).fileNameList = run(iRun).fileNameList(iSet,:);
        end
        
        % and then go on to the next round
        
    end % for iRun = 1:nRuns
    
    % finally, throw away empty newRuns
    newRun(newRunCt+1:end) = [];
    
    % turn on warnings
    warning(oldWarnings);
    
end %  if split/subset


run = newRun;

%===================================================================
%------END COMPARE TWO DATA SETS-----
%===================================================================



%============================================
%                DOWNSAMPLE
%============================================
if downsample - 1 % triggers if downsample is 2 or more
    rmRun = [];
    for iRun = 1:length(run)
        
        data = run(iRun).data;
        fileNameList = run(iRun).fileNameList;
        
        % calculate length of new array
        lengthData    = length(data);
        lengthNewData = lengthData * downsample;
        
        % assign new output. First, get all the fieldnames of data. Then init
        % new structure and new fileNameList
        dataFieldNames = fieldnames(data);
        newData(1:lengthNewData,1) = cell2struct(cell(size(dataFieldNames)),dataFieldNames,1);
        newFileNameList = cell(lengthNewData,1);
        
        newDataCt = 1;
        
        % loop through data. downsample by taking timepoints 1:ds:nTp
        for iData = 1:lengthData
            
            % set range for iNewData. Top row: index into newData. BottomRow:
            % nth downsampling turn
            %iNewDataRange = [1:downsample; ((iData-1)*downsample+1) : iData*downsample]; 
            
            % get the max timepoint of data (data.timePoints is sorted!)
            dataTpMax = data(iData).timePoints(end);
            
            for iNewData = 1:downsample %iNewDataRange
                
                % read the timepoint-vector
                new2oldTimePoints = [iNewData(1):downsample:dataTpMax]';
                
                % find the where the newTimePoints are in the old timePoints
                % vector (if they are there at all!)
                [dummy,newTp2oldTpIdx] = ismember(new2oldTimePoints,data(iData).timePoints);
                
                % get rid of zeros and create new timepoint vector (with the eventual missing indices)
                [newTimePoints, dummy, newTp2oldTpIdx] = find(newTp2oldTpIdx);
                
                if length(newTimePoints) < constants.MINTRAJECTORYLENGTH
                    % we do not save anything
                else
                    
                    newData(newDataCt).timePoints = newTimePoints;% newData(iNewData(2)).timePoints = newTimePoints;
                    
                    % read time, distance
                    newData(newDataCt).distance = data(iData).distance(newTp2oldTpIdx,:);
                    newData(newDataCt).time     = data(iData).time(newTp2oldTpIdx,:);
                    
                    % don't forget info
                    newData(newDataCt).info = data(iData).info;
                    
                    % update fileNameList
                    newFileNameList{newDataCt,1} = [fileNameList{iData} '_ds' num2str(downsample) '_' num2str(iNewData(1))];
                    
                    newDataCt = newDataCt + 1;
                    
                end % if length(newTimePoints) < constants.MINTRAJECTORYLENGTH
                
            end % for iNewData = 1:downsample
            
        end % for iData
        
        % remove empty entries into data/filenamelist
        if newDataCt <= lengthNewData
            newData(newDataCt:end) = [];
            newFileNameList(newDataCt:end) = [];
            if newDataCt == 1 % if all removed
                rmRun = [rmRun;iRun];
            end
        end
        
        % assign output
        run(iRun).data = newData;
        run(iRun).fileNameList = newFileNameList;
        clear newData newFileNameList iData iNewData lengthData lengthNewData dataFieldNames newDatact
    end % for iRun
    
    % remove runs with no data
    run(rmRun) = [];
    if isempty(run)
        error('after downsampling, no data was left at all')
    end
end % if downsample

%============================================
%        END DOWNSAMPLE
%============================================



%========================================
%        RANDOMIZE
%========================================
if randomize
    
    for iRun = 1:length(run)
        
        
        % put data into random order
        nData = length(run(iRun).data);
        randIdx = randperm(nData);
        run(iRun).data = run(iRun).data(randIdx);
        run(iRun).fileNameList = run(iRun).fileNameList(randIdx);
        
    end % for iRun = 1:length(run)
    
end