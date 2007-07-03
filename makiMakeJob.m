function job = makiMakeJob(jobType)
%MAKIMAKEJOB is a hack to set up jobs from movies
%
% SYNOPSIS: makiMakeJob
%
% INPUT
%
% OUTPUT job: job-struct (input to makiMovieAnalysis)
%
% REMARKS This is a hack. Setting up jobs should be done via a GUI.
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

% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER


switch jobType
    case 1
        % test jobs
        jobPath = 'D:\makiTestData';
        

        % find all .dv files in jobPath
        fileList = searchFiles('dv$','log',jobPath,1);

        selectIdx = listSelectGUI(fileList(:,1),[],'move');
        % shorten fileList
        fileList = fileList(selectIdx,:);

        nJobs = length(selectIdx);

        % set up job file
        job(1:nJobs) = struct('jobPath','','jobName','','dataStruct',[]);

        job(1).jobPath = jobPath;
        job(1).jobName = sprintf('testJob-%s.mat',nowString);
        %!!!! for future jobs: also put winRoot, linuxRoot
        % make as cell, so that it is easily accessed via ispc

        for iJob = 1:nJobs
            % for now: project name is rawMovieName minus extension
            rawMovieName = fileList{iJob,1};
            rawMoviePath = fileList{iJob,2};
            [checkPath,projectName] = fileparts(rawMovieName);
            
            % movie could be in a subdir of the official job path
            if strmatch(jobPath,checkPath)
                jobPath = checkPath;
            end

            % make folder for file and movie if necessary
            dataFilePath = fullfile(jobPath,projectName);
            if ~isdir(dataFilePath)
                mkdir(dataFilePath);
                movefile(fullfile(rawMoviePath, rawMovieName),dataFilePath);
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
                dataStruct.dataProperties = defaultDataProperties(dataStruct.dataProperties);
            else
                %bug. do dir.name dataStruct = makiLoadDataFile(dataFile(end,:));
            end

            % crop if not cropped yet (do it here so that we can easily update
            % dataProperties
            if dataStruct.status(1) < 1 && ispc
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
                    'that there are no kinetochores from other cells in the image.',...
                    'Careful: Crop BEFORE clicking OK!'],'Cropping');
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
                imarisHandle = [];
            end


            % request initCoord if necessary
            toDo = find(dataStruct.status(3)<1);
            dataStruct.status(toDo+2) = -1;

            % store dataStruct - store success
            dataStruct.status(2) = makiSaveDataFile(dataStruct);

            % update job
            job(iJob).dataStruct = dataStruct;

        end % loop jobs
    otherwise
        error('job type %i not defined', jobType)
end

% save job
save(fullfile(job(1).jobPath,job(1).jobName),'job');

warning(warningState)
