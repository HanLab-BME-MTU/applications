function job = makiMakeJob(jobType,status,job,ask4input)
%MAKIMAKEJOB is a hack to set up jobs from movies
%
% SYNOPSIS: job = makiMakeJob(jobType,status,job)
%
% INPUT jobType: string which can take the values:
%                'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or 
%                'MCAINSH'
%       status : status that you want to achieve, either as status vector,
%                e.g. [1,1,1,0,1,0,0,0,0], or a list of tasks, e.g. [1,2,3,5]
%                To force a task, set 2 in the status vector, or negative
%                numbers in the list
%               
%                REMARK: this help page (and the code) freely mixes the
%                terms 'job' and 'task'. A job contains a list of tasks
%                that is applied to one movie. Thus, status effectively defines 
%                a task list 
%
%       job    : job-struct as output from makiMakeJob (e.g. if you want to
%                update it). If empty, you can load via GUI
%
%       ask4input: 1 to ask user for input (specifically for tracking and
%                  sister identification), 0 to use built-in defaults.
%                  Optional. Default: 1.
%
% OUTPUT job: job-struct (input to makiMovieAnalysis)
%
% REMARKS Order of tasks so far: 
%         (1) cropping
%         (2) saving dataStruct
%         (3) initCoord
%         (4) mixture-model fitting - DON'T USE!
%         (5) plane fit
%         (6) tracking 
%         (7) sister identification
%         (8) update kinetochore and frame classification
%         (9) frame alignment
%
%         I removed the option jobType = 'UPDATE' for now because of
%         changes that I had to make to makiPathDef due to more than
%         one server being called O: on the Windows side.
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn, kjaqaman, gdanuser
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make this compatible for future addition of tasks.
makiNumPossibleTasks = 9;

if nargin < 1 || isempty(jobType)
    jobType = 1;
    serverType = 'TEST';
else
    jobType = upper(jobType);
    serverType = jobType;
    if strcmp(jobType,'TEST')
        jobType = 1;
    elseif strcmp(jobType,'HERCULES')
        jobType = 2;
    elseif strcmp(jobType,'UPDATE')
        jobType = 3;
    elseif strcmp(jobType,'DANUSER')
        jobType = 4;
    elseif strcmp(jobType,'MERALDI')
        jobType = 5;
    elseif strcmp(jobType,'SWEDLOW')
        jobType = 6;
    elseif strcmp(jobType,'MCAINSH')
        jobType = 7;
    end
end
if nargin < 2 || isempty(status)
    status = zeros(makiNumPossibleTasks,1);
    status(3) = 1; % do only initCoord for now as standard
else
    status = status(:);
end

% Handle the status input as a task list instead of a status vector
if length(status) < makiNumPossibleTasks || any(status > 2)
    tmp = zeros(makiNumPossibleTasks,1);
    tmp(abs(status)) = 1;
    tmp(abs(status(status<0))) = 2;
    status = tmp;
end

if jobType == 3 && nargin < 3
    job = [];
end

if nargin < 4 || isempty(ask4input)
    ask4input = 1;
end

% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER

switch jobType
    case {1,2,4,5,6,7}
        % test jobs and hercules runs

        % basePath depends on the job
        switch jobType
            case 1
                jobPath = [makiPathDef('$TESTDATA',serverType) filesep];
                basePath = jobPath;
                jobName = sprintf('testJob-%s.mat',nowString);
            case 2
                jobPath = [makiPathDef('$TESTDATA',serverType) filesep];
                basePath = [makiPathDef('$HERCULES',serverType) filesep];
                jobName = sprintf('herculesJob-%s.mat',nowString);
            case 4
                basePath = [makiPathDef('$DANUSER',serverType) filesep];
                jobPath = [basePath 'JobLogFiles'];
                jobName = sprintf('danuserJob-%s.mat',nowString);
            case 5
                basePath = [makiPathDef('$MERALDI',serverType) filesep];
                jobPath = [basePath 'JobLogFiles'];
                jobName = sprintf('meraldiJob-%s.mat',nowString);
            case 6
                basePath = [makiPathDef('$SWEDLOW',serverType) filesep];
                jobPath = uigetdir(basePath,'Please select JobLogFiles directory');
                jobName = sprintf('swedlowJob-%s.mat',nowString);
            case 7
                basePath = [makiPathDef('$MCAINSH',serverType) filesep];
                jobPath = [basePath 'JobLogFiles'];
                jobName = sprintf('mcainshJob-%s.mat',nowString);
        end

        % allow user to change base path
        basePath = uigetdir(basePath,'Please select directory of movies to be analyzed');

        % find all .dv files in basePath
        fileList = searchFiles('dv$','(log)|(_PRJ)',basePath,1);

        selectIdx = listSelectGUI(fileList(:,1),[],'move');
        % shorten fileList
        fileList = fileList(selectIdx,:);

        nJobs = length(selectIdx);

        % set up job file
        % store job path for all jobs, in case we want to pick out one
        job(1:nJobs) = struct('jobPath',jobPath,...
            'jobName',jobName,'dataStruct',[]);

        %initialize progress display
        progressText(0,'Setting up jobs');

        for iJob = 1:nJobs

            % for now: project name is rawMovieName minus extension
            rawMovieName = fileList{iJob,1};
            rawMoviePath = fileList{iJob,2};

            % for the project name: remove _R3D
            % works only when _R3D is immediately before the extension
            % At this point (July-24-2007), we have a number of files (and
            % directories) which have the _R3D string in the middle of the
            % name and thus will contain it also in the project name
            extIdx = regexp(rawMovieName,'(_R3D)?\.dv');
            projectName = rawMovieName(1:extIdx-1);

            % make folder for file and movie if necessary
            % For both test and hercules projects, rawMoviePath and
            % dataFilePath are identical. Thus, check whether the
            % rawMoviePath already contains the project name. If not, make
            % a subdirectory
            if any(findstr(rawMoviePath,projectName))
                dataFilePath = rawMoviePath;
            else
                % movie isn't in its directory yet. Make a new one!
                dataFilePath = fullfile(rawMoviePath,projectName);
                mkdir(dataFilePath);
                % move everything that contains the project name (such as
                % the log, or a projection)
                % CORRECTION HERE: Since the _R3D and the crop identifier
                % strings can occur in random positions, the project name
                % sometimes is and sometimes is not the general file name body
                % for the '.dv', '_PRJ.dv', and '.dv.log' files. Thus, the previous call  
                % files2move =
                % searchFiles(projectName,'(log)|(_PRJ)',rawMoviePath,0)
                % resulted in either the moving of all, or of only the
                % movie file. 
                % For now, I restrict the moving to the movie file only,
                % the _PRJ remains untouched, and the '.dv.log' file is
                % dealt with specifically below. 
                rawMovieFile = searchFiles(projectName,'(log)|(_PRJ)',rawMoviePath,0);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD VERSION
                %                 for iFile = 1:size(files2move,1)
                %                     movefile(...
                %                         fullfile(files2move{iFile,2},files2move{iFile,1}),...
                %                         dataFilePath);
                %                 end
                if size(rawMovieFile,1) > 1
                    errordlg(sprintf('%s is not a unique movie filename in the directory %s\nRename before continuing', ...
                        rawMovieFile{1,1}, rawMovieFile{1,2}, ...
                        'Non unique movie filename'));
                    return;
                else
                    movefile(...
                        fullfile(rawMovieFile{1,2},rawMovieFile{1,1}),...
                        dataFilePath)
                end
                
                % NOW: copy log file

                % get 'body' of movie name, e.g. from 'bla_crop1_xx.dv'
                % extract 'bla'; or from 'bla_R3D_crop1.dv' extra 'bla'
                cropIdx = regexp(projectName,'(_CPY)|(_cpy)|(_crop)|(_CROP)');
                if ~isempty(cropIdx)
                    % check if the project name still contains a _R3D
                    r3dIdx = regexp(projectName,'(_R3D)');
                    if ~isempty(r3dIdx)
                        % the body of the project name is from 1 to 
                        % either a '_R3D' or one of the crop strings
                        logFileNameBody = projectName(1:min(cropIdx,r3dIdx));
                    else
                        % the body of the project name is from 1 to 
                        % one of the crop strings
                        logFileNameBody = projectName(1:cropIdx);
                    end
                else
                    logFileNameBody = projectName;
                end;
                % search for bla_*log
                logFile = searchFiles([logFileNameBody,'.*','\.log$'],'',rawMoviePath,0);
                % copy file log to dataFilePath
                if size(logFile,1) > 1
                    errordlg(sprintf('%s is not a unique log filename in the directory %s\nRename before continuing', ...
                        logFile{1,1}, logFile{1,2}, ...
                        'Non unique log filename'));
                else
                    if ~isempty(logFile)
                        copyfile(...
                            fullfile(logFile{1,2},logFile{1,1}),...
                            dataFilePath)
                        % rename the log file in dataFilePath to
                        % movieFileName.log 
                        % there are a few rare cases where the moved log
                        % file already has the movieFileName 
                        % This happens when the movie itself has no crop
                        % indicator
                        movedLogFileOldName = fullfile(dataFilePath,logFile{1,1});
                        movedLogFileNewName = fullfile(dataFilePath,[rawMovieFile{1,1},'.log']);
                        if ~strmatch(movedLogFileOldName,movedLogFileNewName)
                            movefile(movedLogFileOldName,movedLogFileNewName);
                        end
                    else
                        warndlg(sprintf('No log filename matching movie %s in the directory %s\n', ...
                            rawMovieFile{1,1}, rawMovieFile{1,2}), ...
                            'Non unique log filename');
                    end
                end
                
                rawMoviePath = dataFilePath;
            end

            % make dataStruct if necessary
            dataFile = dir(fullfile(dataFilePath,'*-makiData-*'));
            if isempty(dataFile)
                dataStruct = makiMakeDataStruct(rawMovieName, ...
                    rawMoviePath, projectName, dataFilePath);

                % make movieHeader
                dataStruct.movieHeader = readr3dheader(fullfile(rawMoviePath,rawMovieName));

                % make dataProperties, set sigmaCorrection to 1.5
                dataStruct.dataProperties = defaultDataProperties(dataStruct.movieHeader);
                dataStruct.dataProperties.sigmaCorrection = [1.5,1.5];
                dataStruct.dataProperties.MAXSPOTS = 500; % number of locmax considered by initCoord

                % run defDP again to get the correct filter parameters
                dataStruct.dataProperties = defaultDataProperties(dataStruct.dataProperties);

                % set crop-status to -1 if we're on a PC
                if ispc
                    dataStruct.status(1) = -1;
                else
                    dataStruct.status(1) = 0;
                end
            else
                % load dataFile
                dataStruct = makiLoadDataFile(serverType,fullfile(dataFilePath,dataFile.name));
            end

            % request jobs to be done
            toDo = dataStruct.status<status;
            dataStruct.status(toDo) = -1;

            % crop if not cropped yet (do it here so that we can easily update
            % dataProperties, and that we can set up the job on windows and
            % then run under linux)
            if dataStruct.status(1) < 0 && ispc
                % load movie into Imaris
                [dummy,dummy,dummy,...
                    dummy,dummy,imarisHandle] = ...
                    imarisImread(...
                    rawMovieName,...
                    rawMoviePath,[],[],1);

                % read 'zero'. In DV files (and maybe others), Imaris puts
                % the zero at -0.5 pix (like Matlab!). Furthermore, if the
                % movie was cropped already, zero is in the coordinates of
                % the original movie. Therefore, subtract zeroOffset below
                zeroOffsetX = imarisHandle.mDataSet.mExtendMinX;
                zeroOffsetY = imarisHandle.mDataSet.mExtendMinY;
                %zeroOffsetZ = imarisHandle.mDataSet.mExtendMinZ;

                % in principle, we could also allow cropping in Z. In
                % practice, there is no need for that
                announcement = warndlg(['Please crop this movie in xy so',...
                    ' that there are no kinetochores from other cells in the image.',...
                    'CAREFUL: Crop and WAIT till Imaris is done cropping BEFORE clicking OK!'],'Please read carefully');
                uiwait(announcement);

                %%%%%%%%%%%%%%%%% THIS CROPPING WILL BE MODIFIED TO
                %%%%%%%%%%%%%%%%% INTEGRATE TIME CROPPING AS WELL -- gd
                
                % crop has to be corrected by zeroOffset. Round to avoid
                % .9999999 pixels
                % Also, correct for the fact that the image is 0.5 pixel
                % 'too large'!
                % !!!! SWAP X/Y !!!!
                cropYmin = round((imarisHandle.mDataSet.mExtendMinX-zeroOffsetX)/dataStruct.dataProperties.PIXELSIZE_XY)+1;
                cropYmax = round((imarisHandle.mDataSet.mExtendMaxX-zeroOffsetX)/dataStruct.dataProperties.PIXELSIZE_XY)+1-1;
                cropXmin = round((imarisHandle.mDataSet.mExtendMinY-zeroOffsetY)/dataStruct.dataProperties.PIXELSIZE_XY)+1;
                cropXmax = round((imarisHandle.mDataSet.mExtendMaxY-zeroOffsetY)/dataStruct.dataProperties.PIXELSIZE_XY)+1-1;
                %cropZmin = round((imarisHandle.mDataSet.mExtendMinZ-zeroOffsetZ)/dataStruct.dataProperties.PIXELSIZE_Z)+1;
                %cropZmax = round((imarisHandle.mDataSet.mExtendMaxZ-zeroOffsetZ)/dataStruct.dataProperties.PIXELSIZE_Z)+1-1;
                % check whether we're cropping at all

                % GAUDENZ:  check for negative

                if cropXmin == 1 && cropYmin == 1 && ...
                        cropXmax == dataStruct.dataProperties.movieSize(1) &&  ...
                        cropYmax == dataStruct.dataProperties.movieSize(2)
                    % do nothing (or 'update' dataProperties)
                    dataStruct.dataProperties.crop = [];
                else
                    % store crop, adjust movieSize
                    dataStruct.dataProperties.crop = ...
                        [[cropXmin,cropYmin;cropXmax,cropYmax],zeros(2,3)];
                    dataStruct.dataProperties.movieSize(1) = cropXmax - cropXmin + 1;
                    dataStruct.dataProperties.movieSize(2) = cropYmax - cropYmin + 1;
                end

                % set crop status
                dataStruct.status(1) = 1;
                dataStruct.statusHelp{1,2} = date;

                % close imaris to reduce memory consumption (otherwise, two
                % imaris sessions will open at the same time)
                clear imarisHandle
            end

            %assign tracking parameters
            if dataStruct.status(6) < 0
                dataStruct = getTrackingInput(dataStruct,ask4input);
            end

            %assign sister grouping parameters
            if dataStruct.status(7) < 0
                dataStruct = getGroupSisterInput(dataStruct,ask4input);
            end

            % store dataStruct - store success
            [dataStruct.status(2),dataStruct] = makiSaveDataFile(serverType,dataStruct,'dataProperties');
            % dataStruct.status(2) = makiSaveDataFile(dataStruct,'dataProperties');
            % update job
            job(iJob).dataStruct = dataStruct;

            %display progress
            progressText(iJob/(nJobs),'Setting up jobs');

        end % loop jobs
    case 3
        % job file input

        % check input
        if isempty(job)
            [fname,pname] = uigetfile('*job*','select job file');
            % load from file
            if fname == 0
                return
            else
                load(fullfile(pname,fname));
            end
        else
            fname = [];
        end

        % loop through job to reload dataStruct, and to change status
        job = makiMakeJobPlatformIndependent(job,serverType);

        for iJob = 1:length(job)
            % reload dataStruct
            dataStruct = makiLoadDataFile(serverType,...
                fullfile(job(iJob).dataStruct.dataFilePath,...
                job(iJob).dataStruct.dataFileName));

            % set status
            toDo = dataStruct.status<status;
            dataStruct.status(toDo) = -1;

            % ask for input for groupSisters if needed
            if dataStruct.status(7) < 0
                dataStruct = getGroupSisterInput(dataStruct);
            end

            % update job
            job(iJob).dataStruct = dataStruct;
            if ~isempty(fname)
                job(iJob).jobName = fname;
                job(iJob).jobPath = pname;
            end

        end % loop jobs

    otherwise
        error('job type %i not defined', jobType)
end

if ~isempty(job)
    % make job platform-independent
    job = makiMakeJobPlatformIndependent(job,serverType);

    % save job
    save(fullfile(makiPathDef(job(1).jobPath,serverType),job(1).jobName),'job');
end

warning(warningState)



%% subfunction to ask for input for tracking
function dataStruct = getTrackingInput(dataStruct,ask4input)

%define default values
rotate_def = 1;
timeWindow_def = 3;
minRadius_def = 0.5;
maxRadius_def = 1.1;

%assign defaults
dataPropertiesTmp = dataStruct.dataProperties;
if isfield(dataPropertiesTmp,'tracksParam') %if tracksParam have been assigned previously
    tracksParamTmp = dataPropertiesTmp.tracksParam;
    if isfield(tracksParamTmp,'rotate')
        rotateTmp = tracksParamTmp.rotate;
    else
        rotateTmp = rotate_def;
    end
    if isfield(tracksParamTmp,'gapCloseParam')
        timeWindowTmp = tracksParamTmp.gapCloseParam.timeWindow - 1;
    else
        timeWindowTmp = timeWindow_def;
    end
    if isfield(tracksParamTmp,'costMatParam')
        minRadiusTmp = tracksParamTmp.costMatParam.minSearchRadiusL;
        maxRadiusTmp = tracksParamTmp.costMatParam.maxSearchRadiusL;
    else
        minRadiusTmp = minRadius_def;
        maxRadiusTmp = maxRadius_def;
    end
else %if this is the first time
    rotateTmp = rotate_def;
    timeWindowTmp = timeWindow_def;
    minRadiusTmp = minRadius_def;
    maxRadiusTmp = maxRadius_def;
end

%ask for user input
if ask4input
    
    tracksParamIn = inputdlg(...
        {'Use rotated coordinates (1 yes, 0 no)',...
        'Maximum gap to close (in frames)',...
        'Minimum allowed search radius (in microns)',...
        'Maximum allowed search radius (in microns)'},...
        sprintf(['Tracking parameters for ' dataStruct.projectName]),1,...
        {num2str(rotateTmp),num2str(timeWindowTmp),...
        num2str(minRadiusTmp),num2str(maxRadiusTmp)},'on');
    
    if isempty(tracksParamIn)
        error('input aborted')
    else
        rotateTmp = str2double(tracksParamIn{1});
        timeWindowTmp = str2double(tracksParamIn{2});
        minRadiusTmp = str2double(tracksParamIn{3});
        maxRadiusTmp = str2double(tracksParamIn{4});
    end
    
else
    
    rotateTmp = rotate_def;
    timeWindowTmp = timeWindow_def;
    minRadiusTmp = minRadius_def;
    maxRadiusTmp = maxRadius_def;
    
end

%assign whether to use rotated coordinates or not
rotate = rotateTmp;

%assign gap closing parameters
gapCloseParam.timeWindow = timeWindowTmp + 1;
gapCloseParam.mergeSplit = 0;
gapCloseParam.minTrackLen = 1;

%assign cost matrix parameters for linking spots between consecutive 
%frames
costMatParam.minSearchRadiusL = minRadiusTmp;
costMatParam.maxSearchRadiusL = maxRadiusTmp;
costMatParam.brownStdMultL = 3.5;
costMatParam.closestDistScaleL = 2;
costMatParam.maxStdMultL = 100;

%assign cost matrix parameters for closing gaps and (in principle) 
%merging and splitting
costMatParam.minSearchRadiusCG = minRadiusTmp;
costMatParam.maxSearchRadiusCG = maxRadiusTmp;
costMatParam.brownStdMultCG = 3.5*ones(gapCloseParam.timeWindow,1);
costMatParam.linStdMultCG = 3.5*ones(gapCloseParam.timeWindow,1);
costMatParam.timeReachConfB = min(2,gapCloseParam.timeWindow);
costMatParam.timeReachConfL = 1;
costMatParam.closestDistScaleCG = 2;
costMatParam.maxStdMultCG = 100;
costMatParam.lenForClassify = 10;
costMatParam.maxAngleVV = 20;
costMatParam.maxAngleVD = 20;
costMatParam.ampRatioLimitCG = [0.65 4];

%assign parameters for using local density to expand search radius
useLocalDensity.link = 1;
useLocalDensity.cg = 1;
useLocalDensity.nnWindowL = gapCloseParam.timeWindow;
useLocalDensity.nnWindowCG = gapCloseParam.timeWindow;

%save tracking parameters in dataStruct
tracksParam.rotate = rotate;
tracksParam.gapCloseParam = gapCloseParam;
tracksParam.costMatParam = costMatParam;
tracksParam.useLocalDensity = useLocalDensity;
dataStruct.dataProperties.tracksParam = tracksParam;

%% subfunction to ask for input for groupSisters
function dataStruct = getGroupSisterInput(dataStruct,ask4input)

%define default values
useAlignment_def = 1;
maxAngle_def = 30;
maxDist_def = 1.5;
minOverlap_def = 10;
useAnaphase_def = 1;
robust_def = 0;

%assign defaults
dataPropertiesTmp = dataStruct.dataProperties;
if isfield(dataPropertiesTmp,'groupSisters') %if groupSisters parameters have been assigned previously
    groupSistersTmp = dataPropertiesTmp.tracksParam;
    if isfield(groupSistersTmp,'useAlignment')
        useAlignmentTmp = groupSistersTmp.costFuntion;
        if ~isnumeric(useAlignmentTmp)
            useAlignmentTmp = useAlignment_def;
        end
    else
        useAlignmentTmp = useAlignment_def;
    end
    if isfield(groupSistersTmp,'maxAngle')
        maxAngleTmp = groupSistersTmp.maxAngle;
    else
        maxAngleTmp = maxAngle_def;
    end
    if isfield(groupSistersTmp,'maxDist')
        maxDistTmp = groupSistersTmp.maxDist;
    else
        maxDistTmp = maxDist_def;
    end
    if isfield(groupSistersTmp,'minOverlap')
        minOverlapTmp = groupSistersTmp.minOverlap;
    else
        minOverlapTmp = minOverlap_def;
    end
    if isfield(groupSistersTmp,'useAnaphase')
        useAnaphaseTmp = groupSistersTmp.useAnaphase;
    else
        useAnaphaseTmp = useAnaphase_def;
    end
    if isfield(groupSistersTmp,'robust')
        robustTmp = groupSistersTmp.robust;
    else
        robustTmp = robust_def;
    end
else %if this is the first time
    useAlignmentTmp = useAlignment_def;
    maxAngleTmp = maxAngle_def;
    maxDistTmp = maxDist_def;
    minOverlapTmp = minOverlap_def;
    useAnaphaseTmp = useAnaphase_def;
    robustTmp = robust_def;
end

%ask for user input
if ask4input

    groupStats = inputdlg(...
        {'Use alignment (1 yes, 0 no)',...
        'Maximum average angle between sisters and normal to metaphase plate (in degrees)',...
        'Maximum average distance between sisters (in microns)',...
        'Minimum overlap between two tracks (>=10 frames)',...
        'Assume sisters will go in opposite directions in anaphase (1 yes, 0 no)',...
        'Use robust statistics'},sprintf(['Sister identification parameters for ' ...
        dataStruct.projectName]),1,{num2str(useAlignmentTmp),...
        num2str(maxAngleTmp),num2str(maxDistTmp),num2str(minOverlapTmp),...
        num2str(useAnaphaseTmp),num2str(robustTmp)},'on');

    if isempty(groupStats)
        error('input aborted')
    else
        useAlignmentTmp = str2double(groupStats{1});
        maxAngleTmp = str2double(groupStats{2});
        maxDistTmp = str2double(groupStats{3});
        minOverlapTmp = str2double(groupStats{4});
        useAnaphaseTmp = str2double(groupStats{5});
        robustTmp = str2double(groupStats{6});
    end

else

    useAlignmentTmp = useAlignment_def;
    maxAngleTmp = maxAngle_def;
    maxDistTmp = maxDist_def;
    minOverlapTmp = minOverlap_def;
    useAnaphaseTmp = useAnaphase_def;
    robustTmp = robust_def;

end

%process user input
groupSisters.useAlignment = useAlignmentTmp;
groupSisters.maxAngle = maxAngleTmp;
groupSisters.maxDist = maxDistTmp;
groupSisters.minOverlap = minOverlapTmp;
groupSisters.useAnaphase = useAnaphaseTmp;
groupSisters.robust = robustTmp;

%save in dataStruct
dataStruct.dataProperties.groupSisters = groupSisters;
