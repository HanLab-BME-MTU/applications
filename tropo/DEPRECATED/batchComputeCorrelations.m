function batchComputeCorrelations(varargin)
%      BATCHCOMPUTECORRELATIONS(forceRedo) Compute the correlation
%      between:
%      - average Actin speed along Actin tracks (STEPS 2-10)
%      - average TM speed along Actin tracks (STEPS 2-10)
%      - average protrusion value along Actin tracks (STEP 2-10)
%      - cell activity and speckle distance to the edge (STEP 11)
% 'forceRun' forces to recompute results.

if nargin >= 1 && ~isempty(varargin{1})
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if nargin >= 2 && ~isempty(varargin{2})
    forceRun = varargin{2};
else
    forceRun = zeros(11, 1);
end

if nargin >= 3 && ~isempty(varargin{3})
    batchMode = varargin{3};
else
    batchMode = 1;
end

% Get every path from rootDirectory containing two subfolders among the
% specified list of directories.
subDirNames = {'actin', 'TM2', 'TM4', 'TM5NM1'};
paths = getDirectories(rootDirectory, subDirNames, ...
    @(x) exist([x filesep 'lastProjSettings.mat'], 'file'));

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories (Grab a coffee)...');

nMovies = numel(paths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    path = paths{iMovie};
   
    currMovie = movieData{iMovie};
    
    %
    % STEP 1: Create the initial movie data
    %
    
    currMovie.analysisDirectory = [paths{iMovie} filesep 'windowAnalysis'];
    currMovie.channels(1:2) = struct('fsmDirectory', [], 'imageDirectory', []);
    
    % STEP 1.1: Get the FSM directories
    ind = 0;
    for i = 1:numel(subDirNames)
        subDirPath = [path filesep subDirNames{i}];
        if exist(subDirPath, 'dir')
            ind = ind + 1;
            currMovie.channels(ind).fsmDirectory = subDirPath;
        end
    end
    
    if ind ~= 2
        disp([movieName ': Unable to find the 2 FSM directories. (SKIPPING)']);
        continue;
    end

    currMovie.channels(1).imageDirectory = [currMovie.channels(1).fsmDirectory filesep 'crop'];
    currMovie.channels(2).imageDirectory = [currMovie.channels(2).fsmDirectory filesep 'crop'];
    
    nImages1 = numel(dir([currMovie.channels(1).imageDirectory filesep '*.tif']));
    nImages2 = numel(dir([currMovie.channels(2).imageDirectory filesep '*.tif']));
    
    if nImages1 ~= nImages2
        disp([movieName ': Number of images in two channels differ. (SKIPPING)']);
        continue;
    end
    
    currMovie.nImages = nImages1;

    % STEP 1.2: Load physical parameter from
    load([currMovie.channels(1).fsmDirectory filesep 'fsmPhysiParam.mat']);
    currMovie.pixelSize_nm = fsmPhysiParam.pixelSize;
    currMovie.timeInterval_s = fsmPhysiParam.frameInterval;
    clear fsmPhysiParam;

    % STEP 1.3: Get the mask directory
    currMovie.masks.directory = [currMovie.channels(1).fsmDirectory filesep 'edge' filesep 'cell_mask'];
    if ~exist(currMovie.masks.directory, 'dir')
        disp([movieName ': unable to find mask directory. (SKIPPING)']);
        continue;
    end
    currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));
    currMovie.imageDirectory = currMovie.channels(1).imageDirectory;
    
    % STEP 1.4: Update from already saved movieData
    if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && ~forceRun(1)
        currMovie = load([currMovie.analysisDirectory filesep 'movieData.mat']);
        currMovie = currMovie.movieData;
    end
    
    %
    % STEP 2: Compute the average distance of each track to the cell edge.
    %

    if ~isfield(currMovie,'meanDistances') || ~isfield(currMovie.meanDistances,'status') || ...
            currMovie.meanDistances.status ~= 1 || forceRun(2)
        try
            % Format:
            % meanDistances(i, 1) = mean distance of the ith track to the cell edge
            % meanDistances(i, 2) = life time of the ith track
            meanDistances = cell(2, 1);

            currMovie.meanDistances.status = 0;
            currMovie.meanDistances.directory = currMovie.analysisDirectory;
            currMovie.meanDistances.filename = 'meanDistances.mat';
            
            % STEP 2.1: Load MPM file
            MPMs = cell(2, 1);
            
            for k = 1:2
                load([currMovie.channels(k).fsmDirectory filesep 'tack' filesep 'mpm.mat']);                        
                MPMs{k} = MPM;
            end

            clear M MPM;
            
            % STEP 2.2: Compute the distance transforms            
            D = cell(numel(currMovie.masks.n), 1);
            
            filenames = dir([currMovie.masks.directory filesep '*.tif']);
            for iFrame = 1:currMovie.masks.n
                BW = imread([currMovie.masks.directory filesep filenames(iFrame).name]);
                D{iFrame} = single(bwdist(max(BW(:)) - BW));
            end
            
            % STEP 2.3: Compute track mean distance to the edge            
            for k = 1:2
                numTracks = 0;
                meanDistances{k} = zeros(numel(find(MPMs{k} == 0)) / 2 + 1, 2);
                for i = 1:size(MPMs{k}, 1)
                    j = 1;                    
                    while j < size(MPMs{k}, 2) - 1
                        if MPMs{k}(i, j:j+1)
                            birth = ceil(j / 2);
                            numTracks = numTracks + 1;
                            d = 0;
                            m = zeros(1, 2);
                            while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
                                iFrame = ceil(j / 2);
                                p = MPMs{k}(i, j:j+1);
                                death = iFrame;
                                d = d + D{iFrame}(p(1), p(2));
                                m = m + p;
                                j = j + 2;
                            end
                            lifeTime = death - birth + 1;
                            meanDistances{k}(numTracks, 1) = d / lifeTime;
                            meanDistances{k}(numTracks, 2) = lifeTime;
                        end
                        
                        j = j + 2;
                    end
                end
        
                meanDistances{k} = meanDistances{k}(1:numTracks, :);
            end
    
            save([currMovie.meanDistances.directory filesep ...
                currMovie.meanDistances.filename], 'meanDistances');
            
            clear D MPMs meanDistances;
            
            currMovie.meanDistances.status = 1;
            updateMovieData(currMovie);
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.meanDistances.error = errMess;
            currMovie.meanDistances.status = 0;
        end
    end
    
    %
    % STEP 3: Compute the average speed along each track.
    %
    if ~isfield(currMovie,'meanSpeeds') || ~isfield(currMovie.meanSpeeds,'status') || ...
            currMovie.meanSpeeds.status ~= 1 || forceRun(3)
        try
            % Format:
            % meanSpeeds(i, 1) = mean Actin speed of the ith track
            % meanSpeeds(i, 2) = mean TM speed of the ith track
            % meanSpeeds(i, 3) = life time of the ith track
            meanSpeeds = cell(2, 1);

            currMovie.meanSpeeds.status = 0;
            currMovie.meanSpeeds.directory = currMovie.analysisDirectory;
            currMovie.meanSpeeds.filename = 'meanSpeeds.mat';
            
            % STEP 3.1: Load MPM file            
            MPMs = cell(2, 1);
            
            for k = 1:2
                load([currMovie.channels(k).fsmDirectory filesep 'tack' filesep 'mpm.mat']);                        
                MPMs{k} = MPM;
            end

            clear M MPM;
            
            % STEP 3.2: Load the speed maps
            S = cell(currMovie.masks.n, 2);
    
            for k = 1:2
                filenames = dir([currMovie.channels(k).fsmDirectory filesep ...
                    'post' filesep 'mat' filesep 'speedMap*.mat']);
                for iFile = 1:numel(filenames)
                    filename = [currMovie.channels(k).fsmDirectory filesep 'post' filesep ...
                        'mat' filesep filenames(iFile).name];
                    load(filename);
                    [dummy1, dummy2, no] = getFilenameBody(filename);
                    no = str2double(no);
                    S{no, k} = single(speedMap);
                    clear speedMap;
                end
                % Extend the first
                filename = [currMovie.channels(k).fsmDirectory filesep ...
                    'post' filesep 'mat' filenames(1).name];
                [dummy1, dummy2, no] = getFilenameBody(filename);
                no = str2double(no);
                for i = 1:no - 1
                    S{i, k} = S{no, k};
                end
                % Extend the last
                filename = [currMovie.channels(k).fsmDirectory filesep ...
                    'post' filesep 'mat' filenames(end).name];
                [dummy1, dummy2, no] = getFilenameBody(filename);
                no = str2double(no);
                for i = no+1:size(S, 1)
                    S{i, k} = S{no, k};
                end
            end

            % STEP 3.3: Compute track mean speed
            for k = 1:2
                numTracks = 0;
                meanSpeeds{k} = zeros(numel(find(MPMs{k} == 0)) / 2 + 1, 3);
                for i = 1:size(MPMs{k}, 1)                  
                    j = 1;                    
                    while j < size(MPMs{k}, 2) - 1
                        if MPMs{k}(i, j:j+1)
                            birth = ceil(j / 2);
                            numTracks = numTracks + 1;
                            s1 = 0; s2 = 0;
                            while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
                                iFrame = ceil(j / 2);
                                p = MPMs{k}(i, j:j+1);
                                death = iFrame;
                                s1 = s1 + S{iFrame, 1}(p(1), p(2));
                                s2 = s2 + S{iFrame, 2}(p(1), p(2));
                                j = j + 2;
                            end
                            lifeTime = death - birth + 1;
                            meanSpeeds{k}(numTracks, 1) = s1 / lifeTime;
                            meanSpeeds{k}(numTracks, 2) = s2 / lifeTime;
                            meanSpeeds{k}(numTracks, 3) = lifeTime;
                        end
                        j = j + 2;
                    end
                end
            end
    
            save([currMovie.meanSpeeds.directory filesep ...
                currMovie.meanSpeeds.filename], 'meanSpeeds');

            clear S MPMs meanSpeeds;
            
            currMovie.meanSpeeds.status = 1;
            updateMovieData(currMovie);
        
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.meanSpeeds.error = errMess;
            currMovie.meanSpeeds.status = 0;
        end
    end

    %
    % STEP 4: Get the contour
    %
    
    dContour = 15; % ~ 1um
    dWin = 10;
    iStart = 2;
    iEnd = 10;
    winMethod = 'e';    
    
    if ~isfield(currMovie,'contours') || ~isfield(currMovie.contours,'status') || ...
            currMovie.contours.status ~= 1 || forceRun(4)
        try
            currMovie = getMovieContours(currMovie, 0:dContour:500, 0, 0, ...
                ['contours_'  num2str(dContour) 'pix.mat'], batchMode);
        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.contours.error = errMess;
            currMovie.contours.status = 0;
        end
    end

    %
    % STEP 5: Calculate protusion
    %
    if ~isfield(currMovie,'protrusion') || ~isfield(currMovie.protrusion,'status') || ...
            currMovie.protrusion.status ~= 1 || forceRun(5)
        try
            currMovie.protrusion.status = 0;

            currMovie = setupMovieData(currMovie);

            handles.batch_processing = batchMode;
            handles.directory_name = [currMovie.masks.directory];
            handles.result_directory_name = [currMovie.masks.directory];
            handles.FileType = '*.tif';
            handles.timevalue = currMovie.timeInterval_s;
            handles.resolutionvalue = currMovie.pixelSize_nm;
            handles.segvalue = 30;
            handles.dl_rate = 30;

            %run it
            [OK,handles] = protrusionAnalysis(handles);

            if ~OK
                currMovie.protrusion.status = 0;
            else
                if isfield(currMovie.protrusion,'error')
                    currMovie.protrusion = rmfield(currMovie.protrusion,'error');
                end
                
                %currMovie.protrusion.directory = [currMovie.masks.directory];
                % Workaround:
                currMovie.protrusion.directory = [currMovie.masks.directory filesep ...
                    'analysis_dl' num2str(handles.dl_rate)];
                
                currMovie.protrusion.fileName = 'protrusion.mat';
                currMovie.protrusion.nfileName = 'normal_matrix.mat';
                currMovie.protrusion.status = 1;
            end
            
            updateMovieData(currMovie);
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.protrusion.error = errMess;
            currMovie.protrusion.status = 0;
        end
    end

    %
    % STEP 6: Create windows
    %
    
    windowString = [num2str(dContour) 'by' num2str(dWin) 'pix_' num2str(iStart) '_' num2str(iEnd)];

    if ~isfield(currMovie,'windows') || ~isfield(currMovie.windows,'status')  || ...
            currMovie.windows.status ~= 1 || forceRun(6)
        try
            currMovie = setupMovieData(currMovie);

            disp(['Windowing movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieWindows(currMovie,winMethod,dWin,[],iStart,iEnd,[],[],...
                ['windows_' winMethod '_' windowString '.mat'], batchMode);

            if isfield(currMovie.windows,'error')
                currMovie.windows = rmfield(currMovie.windows,'error');
            end

        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.windows.error = errMess;
            currMovie.windows.status = 0;
        end
    end

    %
    % STEP 7: Sample the protrusion vector in each window
    %
    if ~isfield(currMovie,'protrusion') || ~isfield(currMovie.protrusion,'samples') || ...
            ~isfield(currMovie.protrusion.samples,'status') || ...
            currMovie.protrusion.samples.status ~= 1 || forceRun(7)
        try
            disp(['Sampling protrusion in movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieProtrusionSamples(currMovie,['protSamples_' ...
                winMethod '_' windowString  '.mat'],10,100, batchMode);
            
            if isfield(currMovie.protrusion.samples,'error')
               currMovie.protrusion.samples = rmfield(currMovie.protrusion.samples,'error');
           end
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);           
            currMovie.protrusion.samples.error = errMess;
            currMovie.protrusion.samples.status = 0;
        end
        
    end 
    
    %
    % STEP 8: Split the windows into different files.
    %
    if ~isfield(currMovie, 'windows') || ~isfield(currMovie.windows, 'splitted') || ...
            currMovie.windows.splitted ~= 1
        splitWindowFrames(currMovie, [currMovie.analysisDirectory filesep 'windows'], batchMode);
        currMovie.windows.splitted = 1;
    end    
    
    %
    % STEP 9: Create the window labels.
    %
    if ~isfield(currMovie, 'labels') || ~isfield(currMovie.labels, 'status') || ...
            currMovie.labels.status ~= 1 || forceRun(9)
        try
            currMovie = setupMovieData(currMovie);

            disp(['Labeling movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            
            currMovie = getMovieLabels(currMovie, batchMode);

            if isfield(currMovie.labels,'error')
                currMovie.labels = rmfield(currMovie.labels,'error');
            end

        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.labels.error = errMess;
            currMovie.labels.status = 0;
        end            
    end
    
    %
    % STEP 10: Compute average protusion score along each track.
    %    
    if ~isfield(currMovie, 'meanProtrusion') || ~isfield(currMovie.meanProtrusion, 'status') || ...
            currMovie.meanProtrusion.status ~= 1 || forceRun(10)
        try
            % Format:
            % meanProtrusions(i, 1) = mean protrusion score of the ith track
            % meanProtrusions(i, 2) = life time of the ith track
            meanProtrusions = cell(2, 1);
            
            currMovie.meanProtrusions.status = 0;
            currMovie.meanProtrusions.directory = currMovie.analysisDirectory;
            currMovie.meanProtrusions.filename = 'meanProtrusions.mat';
            
            % STEP 10.1: Load labels
            L = cell(currMovie.labels.nFrames, 1);
            filenames = dir([currMovie.labels.directory filesep '*.tif']);
            for iFrame = 1:currMovie.labels.nFrames
                L{iFrame} = imread([currMovie.labels.directory filesep filenames(iFrame).name]);
            end
    
            % STEP 10.2: Load protrusion samples file
            load([currMovie.protrusion.directory filesep ...
                currMovie.protrusion.samples.fileName]);
                
            % Add another frame
            avgProt = horzcat(protrusionSamples.averageNormalComponent,...
                protrusionSamples.averageNormalComponent(:, end));
            clear protrusionSamples;
    
            % STEP 10.3: Load MPM file
            MPMs = cell(2, 1);
            
            for k = 1:2
                load([currMovie.channels(k).fsmDirectory filesep 'tack' filesep 'mpm.mat']);                        
                MPMs{k} = MPM;
            end

            clear M MPM;
            
            % STEP 10.4: Compute mean protrusion score along track
            for k = 1:2
                numTracks = 0;
                meanProtrusions{k} = zeros(numel(find(MPMs{k} == 0)) / 2 + 1, 2);
                for i = 1:size(MPMs{k}, 1)                    
                    j = 1;
                    while j < size(MPMs{k}, 2) - 1
                        if MPMs{k}(i, j:j+1)
                            numTracks = numTracks + 1;
                            prot = 0;
                            count = 0;
                            while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
                                iFrame = ceil(j / 2);
                                p = MPMs{k}(i, j:j+1);
                                iSlice = L{iFrame}(p(1), p(2));
                                if iSlice
                                    prot = prot + avgProt(iSlice, iFrame);
                                    count = count + 1;
                                end
                                j = j + 2;
                            end
                            if count
                                meanProtrusions{k}(numTracks, 1) = prot / count;
                            else
                                meanProtrusions{k}(numTracks, 1) = NaN;
                            end
                            meanProtrusions{k}(numTracks, 2) = count;
                        end
                        j = j + 2;
                    end
                end
            end
    
            save([currMovie.meanProtrusions.directory filesep ...
                currMovie.meanProtrusions.filename], 'meanProtrusions');
            
            clear L MPMs meanProtrusions;
            
            currMovie.meanProtrusion.status = 1;
            updateMovieData(currMovie);          
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.meanProtrusion.error = errMess;
            currMovie.meanProtrusion.status = 0;
        end
    end
    
    %
    % STEP 11: Compute correlation between cell activity and speckle
    % distance to cell edge
    %    
    if ~isfield(currMovie, 'activityVSdistance') || ~isfield(currMovie.activityVSdistance, 'status') || ...
            currMovie.activityVSdistance.status ~= 1 || forceRun(11)
        try
            % Format:
            % activityVSdistance(1).name
            % activityVSdistance(1).activity
            % activityVSdistance(1).distance
            activityVSdistance(1:2) = struct('name', [], 'activity', [], 'distance', []);
            
            currMovie.activityVSdistance.status = 0;
            currMovie.activityVSdistance.directory = currMovie.analysisDirectory;
            currMovie.activityVSdistance.filename = 'activityVSdistance.mat';
            
            % STEP 11.1: Load MPM file
            MPMs = cell(2, 1);
            
            for k = 1:2
                load([currMovie.channels(k).fsmDirectory filesep 'tack' filesep 'mpm.mat']);
                MPMs{k} = MPM;
            end
            
            clear M MPM;
            
            % STEP 11.2: Compute the distance transforms
            D = cell(numel(currMovie.masks.n), 1);
            
            filenames = dir([currMovie.masks.directory filesep '*.tif']);
            for iFrame = 1:currMovie.masks.n
                BW = imread([currMovie.masks.directory filesep filenames(iFrame).name]);
                D{iFrame} = single(bwdist(max(BW(:)) - BW));
            end
            
            % STEP 11.3: Load labels
            L = cell(currMovie.labels.nFrames, 1);
            filenames = dir([currMovie.labels.directory filesep '*.tif']);
            for iFrame = 1:currMovie.labels.nFrames
                L{iFrame} = imread([currMovie.labels.directory filesep filenames(iFrame).name]);
            end
            
            % STEP 11.4: Load protrusion samples file
            load([currMovie.protrusion.directory filesep ...
                currMovie.protrusion.samples.fileName]);
            
            % Add another frame
            avgProt = horzcat(protrusionSamples.averageNormalComponent,...
                protrusionSamples.averageNormalComponent(:, end));
            clear protrusionSamples;
            
            % STEP 11.5: Store correlation
            for k = 1:2
                [dummy, name] = getFilenameBody(currMovie.channels(k).fsmDirectory);
                activityVSdistance(k).name = name;
                nSpeckles = numel(find(MPMs{k})) / 2;
                activityVSdistance(k).activity = zeros(nSpeckles, 1);
                activityVSdistance(k).distance = zeros(nSpeckles, 1);
                
                pos = 1;
                
                for iFrame = 1:currMovie.masks.n
                    t = 2 * iFrame - 1;
                    Xt = MPMs{k}(:, t:t+1);
                    % Get rid of 0 rows
                    Xt = Xt(Xt(:, 1) & Xt(:, 2), :);
                    ind = sub2ind(size(D{iFrame}), Xt(:, 1), Xt(:, 2));
                    iSlice = L{iFrame}(ind);
                    % Get rid of index pointed to label == 0
                    ind = ind(iSlice ~= 0);
                    iSlice = L{iFrame}(ind);
                    % Get rid of index pointed to avgProt == 0
                    ind = ind(~isnan(avgProt(iSlice, iFrame)));
                    iSlice = L{iFrame}(ind);
                    activityVSdistance(k).activity(pos:pos+numel(ind)-1) = avgProt(iSlice, iFrame);
                    activityVSdistance(k).distance(pos:pos+numel(ind)-1) = D{iFrame}(ind);
                    pos = pos + numel(ind);
                end
                
                % Remove trailing zeros in activity and distance arrays
                last = find(activityVSdistance(k).activity(end:-1:1) == 0, 1, 'last');
                activityVSdistance(k).activity = ...
                    activityVSdistance(k).activity(1:end-last);
                activityVSdistance(k).distance = ...
                    activityVSdistance(k).distance(1:end-last);
            end
            save([currMovie.activityVSdistance.directory filesep ...
                currMovie.activityVSdistance.filename], 'activityVSdistance');
            
            clear D L MPMs activityVSdistance;
            
            currMovie.activityVSdistance.status = 1;
            updateMovieData(currMovie);
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.activityVSdistance.error = errMess;
            currMovie.activityVSdistance.status = 0;
        end
    end
    
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMov) ': ' errMess.message], mfileName());
    end
    
    movieData{iMovie} = currMovie;
    
    disp([movieName ': DONE']);
end

end