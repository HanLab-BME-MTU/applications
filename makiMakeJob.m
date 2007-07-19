function job = makiMakeJob(jobType,status,job)
%MAKIMAKEJOB is a hack to set up jobs from movies
%
% SYNOPSIS: makiMakeJob
%
% INPUT jobType: 1: test job (default)
%                2: hercules run
%                3: update existing job
%       status : status that you want to achieve, either as status vector,
%                e.g. [1,1,1,0,1,0,0], or a list of jobs, e.g. [1,2,3,5]
%                To force the job, set 2 in the status vector, or negative
%                numbers in the list
%       job    : job-struct as output from makiMakeJob (e.g. if you want to
%                update it). If empty, you can load via GUI
%
% OUTPUT job: job-struct (input to makiMovieAnalysis)
%
% REMARKS This is a hack. Setting up jobs should be done via a GUI.
%         Order of jobs so far: (1) cropping; (2) saving dataStruct;
%         (3) initCoord; (4) plane fit; (5) tracking; (6) sister
%         identification
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(jobType)
    jobType = 1;
end
if nargin < 2 || isempty(status)
    status = [0,0,1,0,0,0,0]'; % do only initCoord for now as standard
else
    status = status(:);
end
if length(status) < 7 || any(status > 2)
    tmp = zeros(7,1);
    tmp(abs(status)) = 1;
    tmp(abs(status(status<0))) = 2;
    status = tmp;
end
if jobType == 3 && nargin < 3
    job = [];
end

% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER

switch jobType
    case {1,2}
        % test jobs and hercules runs

        % job-files will always be stored in testData
        jobPath = makiPathDef('$TESTDATA');

        % basePath depends on the job
        switch jobType
            case 1
                basePath = jobPath;
                jobName = sprintf('testJob-%s.mat',nowString);
            case 2
                basePath = makiPathDef('$HERCULES');
                jobName = sprintf('herculesJob-%s.mat',nowString);
        end

        % allow user to change base path
        basePath = uigetdir(basePath,'Please select data-dir');

        % find all .dv files in jobPath
        fileList = searchFiles('dv$','(log)|(_PRJ)',basePath,1);

        selectIdx = listSelectGUI(fileList(:,1),[],'move');
        % shorten fileList
        fileList = fileList(selectIdx,:);

        nJobs = length(selectIdx);

        % set up job file
        % store job path in all jobs, in case we want to pick out one
        job(1:nJobs) = struct('jobPath',jobPath,...
            'jobName',jobName,'dataStruct',[]);

        for iJob = 1:nJobs
            % for now: project name is rawMovieName minus extension
            rawMovieName = fileList{iJob,1};
            rawMoviePath = fileList{iJob,2};

            % for the project name: remove _R3D
            extIdx = regexp(rawMovieName,'(_R3D)?.dv');
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
                files2move = searchFiles(projectName,'',rawMoviePath,0);
                for iFile = 1:size(files2move,1)
                    movefile(...
                        fullfile(files2move{iFile,2},files2move{iFile,1}),...
                        dataFilePath);
                end
                % copy log file
                % get 'body' of movie name, e.g. from bla_crop1_xx.dv,
                % extract bla
                % search for bla*log
                % copyfile log to dataFilePath
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

                % add groupSisters - properties
                groupSisters.maxDist = 4; % maximum average sister separation in um
                groupSisters.goodTrackRatio = 0.75; % minimum relative track length for grouping
                groupSisters.robust = false; % whether or not to use robust statistics for determining cost function parameters
                groupSisters.costFunction = 'metaphase'; % cost function type.
                dataStruct.dataProperties.groupSisters = groupSisters;

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
                dataStruct = makiLoadDataFile(fullfile(dataFilePath,dataFile.name));
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

            %assign tracking parameters if tracking is requested
            if dataStruct.status(5) < 0

                %ask for some parameters from user
                userEntryR = input('Enter 1 to use rotated coordinates in tracking. ');
                if userEntryR ~= 1
                    userEntryR = 0;
                end

                disp('Please input - as a row vector - values for the maximum gap');
                disp('(in frames), minimum search radius (in microns) and maximum ');
                disp('search radius (in microns).');
                userEntry = input('');

                %gap closing parameters
                gapCloseParam.timeWindow = userEntry(1);
                gapCloseParam.mergeSplit = 0;

                %cost matrix parameters for linking spots between
                %consecutive frames
                costMatParam.minSearchRadiusL = userEntry(2);
                costMatParam.maxSearchRadiusL = userEntry(3);
                costMatParam.brownStdMultL = 3;
                costMatParam.closestDistScaleL = 2;
                costMatParam.maxStdMultL = 20;

                %cost matrix parameters for closing gaps and (in principle)
                %merging and splitting
                costMatParam.minSearchRadiusCG = userEntry(2);
                costMatParam.maxSearchRadiusCG = userEntry(3);
                costMatParam.brownStdMultCG = 3*ones(gapCloseParam.timeWindow,1);
                costMatParam.linStdMultCG = 3*ones(gapCloseParam.timeWindow,1);
                costMatParam.timeReachConfB = gapCloseParam.timeWindow;
                costMatParam.timeReachConfL = 1;
                costMatParam.closestDistScaleCG = 2;
                costMatParam.maxStdMultCG = 20;
                costMatParam.lenForClassify = 10;
                costMatParam.maxAngle = 10;
                costMatParam.ampRatioLimitCG = [0.5000 2];

                %parameters for using local density to expand search radius
                useLocalDensity.link = 1;
                useLocalDensity.cg = 1;
                useLocalDensity.nnWindowL = gapCloseParam.timeWindow;
                useLocalDensity.nnWindowCG = gapCloseParam.timeWindow;

                %save tracking parameters in dataProperties
                dataStruct.dataProperties.tracksParam.rotate = userEntryR;
                dataStruct.dataProperties.tracksParam.gapCloseParam = gapCloseParam;
                dataStruct.dataProperties.tracksParam.costMatParam = costMatParam;
                dataStruct.dataProperties.tracksParam.useLocalDensity = useLocalDensity;

            end %(if dataStruct.status(5) < 0)

            % allow to change settings (yadda yadda GUI yadda yadda)
            if dataStruct.status(6) < 0
                dataStruct = getGroupSisterInput(dataStruct);
            end

            % store dataStruct - store success
            dataStruct.status(2) = makiSaveDataFile(dataStruct);

            % update job
            job(iJob).dataStruct = dataStruct;

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
        job = makiMakeJobPlatformIndependent(job);

        for iJob = 1:length(job)
            % reload dataStruct
            dataStruct = makiLoadDataFile(...
                fullfile(job(iJob).dataStruct.dataFilePath,...
                job(iJob).dataStruct.dataFileName));

            % set status
            toDo = dataStruct.status<status;
            dataStruct.status(toDo) = -1;

            % ask for input for groupSisters if needed
            if dataStruct.status(6) < 0
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
    job = makiMakeJobPlatformIndependent(job);

    % save job
    save(fullfile(makiPathDef(job(1).jobPath),job(1).jobName),'job');
end

warning(warningState)



% subfunction to ask for input for groupSisters
function dataStruct = getGroupSisterInput(dataStruct)

phases ={'prophase';'prometaphase';'metaphase';'anaphase'};
defPhase = strmatch(...
    dataStruct.dataProperties.groupSisters.costFunction,phases);
groupStats = inputdlg(...
    {sprintf(['Main cell cycle phase',...
    '\n1: Prophase\n2: Prometaphase',...
    '\n3: Metaphase\n4: Anaphase']),...
    'distance cutoff (in microns)',...
    'minimum track length relative to movie length',...
    'use robust statistics'},...
    'Grouping Sisters',1,{num2str(defPhase),...
    num2str(dataStruct.dataProperties.groupSisters.maxDist),...
    num2str(dataStruct.dataProperties.groupSisters.goodTrackRatio),...
    num2str(dataStruct.dataProperties.groupSisters.robust)});
if isempty(groupStats)
    error('input aborted')
end

dataStruct.dataProperties.groupSisters.costFunction = ...
    phases{str2double(groupStats{1})};
dataStruct.dataProperties.groupSisters.maxDist = ...
    str2double(groupStats{2});
dataStruct.dataProperties.groupSisters.goodTrackRatio = ...
    str2double(groupStats{3});
dataStruct.dataProperties.groupSisters.robust = ...
    str2double(groupStats{4});