function [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam,saveFile,redoPerMovieAnalysis,parallel)
%RESULTSINDTIMECOURSE compiles the results of a group of movies making one or multiple timecourse datasets and orders them based on time
%
%SYNOPSIS [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam)
%
%INPUT  ML       : MovieList object containing all movies, either all
%                  belonging to one timecourse, or potentially belonging
%                  to multiple timecourses. The number of timecourses is
%                  indicated in the next variable.
%       caseParam: Structure array with number of entries = number of
%                  timecourses among which the movies will be divided. For
%                  each time course, it contains the following fields:
%           .indx    : Indices of movies in ML belonging to this
%                      timecourse.
%           .timeList: The time point of each movie, with order following
%                      the movie order in ML.
%           .name    : Case name, to be used in naming output variables.
%           .indx0min: Index of movie to be considered at time 0, e.g. when
%                      a ligand is added. If > 1, movies before will have a
%                      negative relative time, movies after wil have a
%                      positive relative time. If 1, relative time and
%                      absolute time are the same.
%        saveFile: Logical that determines if this function saves a file.
%                  The default is 'true'. 0 or 1 instead of true or false
%                  will work.
%        redoAnalysis: Logical that determines if analysis should be redone
%                  per movie if prior analysis per movie exists.
%                  The default is true. 0 or 1 instead of true or false
%                  will work.
%        parallel    : char that is either
%                      'none' (default): Run in a single thread
%                      'cluster'       : Run using parallel pool
%                      'slurm'         : Queue on cluster using sbatch
%                      'load_only'     : Load saved data only
%    
%OUTPUT caseTimeList  : 2-column vector indicating the movie times. Column 1
%                       shows absolute time, Column 2 shows relative time.
%       caseResSummary: Structure array storing various results for each
%                       movie, ordered based on their time (corresponding
%                       to caseTimeList). Contains the fields: 
%           .diffSummary         : Diffusion analysis summary, as output
%                                  by summarizeDiffAnRes.
%           .diffCoefMeanPerClass: Mean diffusion coefficient per motion
%                                  class. Order: Immobile, confined, free,
%                                  directed, undetermined.
%           .confRadMeanPerClass : Mean confinement radius per motion
%                                  class. Same order as above.
%           .ampMeanPerClass     : Mean particle amplitude per motion
%                                  class. Rows in same order as above.
%                                  Columns: first for "absolute" amplitude,
%                                  second for amplitude normalized by
%                                  monomer amplitude, as derived from modal
%                                  analysis of particles in last 20 frames
%                                  of each movie.
%           .ampStatsF20         : Amplitude statistics in first 20
%                                  frame of each movie. 
%                                  Order: mean amplitude, first mode
%                                  mean, first mode std, first mode
%                                  fraction, number of modes, mean
%                                  normalized by monomer amplitude (see
%                                  ampMeanPerClass for details).
%           .ampStatsL20         : Same as ampStatsF20, but for last 20
%                                  frames of movie.
%
%Khuloud Jaqaman, March 2015

%% Input and pre-processing

if nargin < 1
    error('resultsIndTimeCourse: Too few input arguments')
end
% For a single parameter, see end of file
if nargin < 2
    caseParam = [];
end
if nargin < 3
    saveFile = true;
end
if nargin < 4
    redoPerMovieAnalysis = true;
end
if nargin < 5
    parallel = 'none';
end

%get number of movies and number of cases
numMovies = length(ML.movieDataFile_);
numCases = length(caseParam);



%reserve memory for individual movie results
resSummary = struct('diffSummary',[],'diffCoefMeanPerClass',[],...
    'confRadMeanPerClass',[],'ampMeanPerClass',[],'ampStatsF20',[],...
    'ampStatsL20',[],'statsMS',[]);
resSummary = repmat(resSummary,numMovies,1);

%define directory for saving results
dir2save = [ML.movieListPath_ filesep 'analysisKJ'];
if(~exist(dir2save,'dir'))
    mkdir(dir2save);
end

%% Calculations

switch(parallel)
    case 'none'
        %go over all movies
        for iM = 1 : numMovies
            file2savePerMovie = fullfile(dir2save,sprintf([dir2save filesep 'resSummary_movie_%03d.mat'],iM));
            if(redoPerMovieAnalysis || ~exist(file2savePerMovie,'file'))
                resSummary(iM) = resultsIndTimeCoursePerMovie(ML.movieDataFile_{iM},file2savePerMovie);
            else
                saved = load(file2savePerMovie);
                resSummary(iM) = saved.resSummary;
            end

            %read in raw results
            %     disp([num2str(iM) '  ' ML.movieDataFile_{iM}]);
        %     MD = MovieData.load(ML.movieDataFile_{iM}); 
        %% mkitti: Moved to resultsIndTimeCoursePerMovie
        %     iProcDiff = MD.getProcessIndex('MotionAnalysisProcess',1,0); %diffusion analysis and tracks
        %     load(MD.processes_{iProcDiff}.outFilePaths_{1});
        %     iProcMask = MD.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
        %     if ~isempty(iProcMask)
        %         mask = imread(fullfile(MD.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
        %     else
        %         mask = [];
        %     end
        %     
        %     %limit analysis to tracks in mask if supplied
        %     if ~isempty(mask) && any(mask(:)==0)
        %         %keep only tracks in mask
        %         numTracks = length(tracks);
        %         keepTrack = ones(numTracks,1);
        %         for iTrack = 1 : numTracks
        %             xCoord = tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
        %             yCoord = tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
        %             meanPosX = round(nanmean(xCoord(:)));
        %             meanPosY = round(nanmean(yCoord(:)));
        %             keepTrack(iTrack) = mask(meanPosY,meanPosX);
        %         end
        %         indxKeep = find(keepTrack);
        %         tracks = tracks(indxKeep);
        %         diffAnalysisRes = diffAnalysisRes(indxKeep);
        %         %redo diffusion analysis summary
        %         if length(indxKeep) < numTracks
        %             minTrackLen = 5;
        %             probDim = 2;
        %             extractType = 1;
        %             [probMotionType,motionChar] = summarizeDiffAnRes(tracks,minTrackLen,probDim,diffAnalysisRes,extractType);
        %             diffAnalysisRes(1).summary.probMotionType = probMotionType;
        %             diffAnalysisRes(1).summary.motionChar = motionChar;
        %         end
        %     end
        %     
        %     %diffusion analysis summary
        %     diffSummary = diffAnalysisRes(1).summary;
        %     
        %     %diffusion analysis per track
        %     trajClass = vertcat(diffAnalysisRes.classification);
        %     trajClass = trajClass(:,2);
        %     trajDiffCoef = catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef');
        %     trajConfRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');
        %     
        %     %amplitude matrix
        %     tracksMat = convStruct2MatIgnoreMS(tracks);
        %     ampMat = tracksMat(:,4:8:end);
        %     
        %     %amplitude per track
        %     ampMeanPerTraj = nanmean(ampMat,2);
        %     
        %     %average properties per motion class
        %     [ampMeanPerClass,diffCoefMeanPerClass,confRadMeanPerClass] = deal(NaN(5,1));
        %     for i = 1 : 2
        %         indxClass = find(trajClass == (i-1));
        %         ampMeanPerClass(i) = mean(ampMeanPerTraj(indxClass));
        %         diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
        %         confRadMeanPerClass(i) = mean(trajConfRad(indxClass));
        %     end
        %     for i = 3 : 4
        %         indxClass = find(trajClass == (i-1));
        %         ampMeanPerClass(i) = mean(ampMeanPerTraj(indxClass));
        %         diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
        %     end
        %     ampMeanPerClass(5) = mean(ampMeanPerTraj(isnan(trajClass)));
        %     
        %     %amplitude statistics in first 20 frames
        %     ampVec = ampMat(:,1:20);
        %     ampVec = ampVec(~isnan(ampVec));
        %     [~,~,modeParam] = fitHistWithGaussians(ampVec,0.05,0,3,0,[],2,[],1,[],0);
        %     numMode = size(modeParam,1);
        %     modeParamMean = modeParam(1,1);
        %     modeParamStd  = modeParam(1,2);
        %     ampMode1Mean = exp( modeParamMean + modeParamStd.^2/2 );
        %     ampMode1Std  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
        %     ampMode1Frac = modeParam(1,4)/sum(modeParam(:,4));
        %     ampStatsF20 = [mean(ampVec) ampMode1Mean ampMode1Std ampMode1Frac numMode];
        %             
        %     %amplitude statistics in last 20 frames
        %     ampVec = ampMat(:,end-19:end);
        %     ampVec = ampVec(~isnan(ampVec));
        %     [~,~,modeParam] = fitHistWithGaussians(ampVec,0.05,0,3,0,[],2,[],1,[],0);
        %     numMode = size(modeParam,1);
        %     modeParamMean = modeParam(1,1);
        %     modeParamStd  = modeParam(1,2);
        %     ampMode1Mean = exp( modeParamMean + modeParamStd.^2/2 );
        %     ampMode1Std  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
        %     ampMode1Frac = modeParam(1,4)/sum(modeParam(:,4));
        %     ampStatsL20 = [mean(ampVec) ampMode1Mean ampMode1Std ampMode1Frac numMode];
        %     
        %     %normalize mean amplitudes by first mode mean from last 20 frames
        %     ampMeanPerClass = [ampMeanPerClass ampMeanPerClass/ampStatsL20(2)]; %#ok<AGROW>
        %     ampStatsF20 = [ampStatsF20 ampStatsF20(1)/ampStatsL20(2)]; %#ok<AGROW>
        %     ampStatsL20 = [ampStatsL20 ampStatsL20(1)/ampStatsL20(2)]; %#ok<AGROW>
        % 
        %     %merge and split statistics
        %     statsMS = calcStatsMS_noMotionInfo(tracks,5,1,1);
        %     
        %     %results for output
        %     resSummary(iM).diffSummary = diffSummary;
        %     resSummary(iM).diffCoefMeanPerClass = diffCoefMeanPerClass;
        %     resSummary(iM).confRadMeanPerClass = confRadMeanPerClass;
        %     resSummary(iM).ampMeanPerClass = ampMeanPerClass;
        %     resSummary(iM).ampStatsF20 = ampStatsF20;
        %     resSummary(iM).ampStatsL20 = ampStatsL20;
        %     resSummary(iM).statsMS = statsMS;

        end
    case 'cluster'
        movieFiles = ML.movieDataFile_;
        parfor iM = 1 : numMovies
            file2savePerMovie = fullfile(dir2save,sprintf([dir2save filesep 'resSummary_movie_%03d.mat'],iM));
            if(redoPerMovieAnalysis || ~exist(file2savePerMovie,'file'))
                resSummary(iM) = resultsIndTimeCoursePerMovie(movieFiles{iM},file2savePerMovie);
            else
                saved = load(file2savePerMovie);
                resSummary(iM) = saved.resSummary;
            end
        end
    case 'slurm'
        scriptName = mfilename('fullpath');
        ML.runSystemCmdPerMovie(['echo ' scriptName 'PerMovie.sh'],@(iM,movieDataFile_) sprintf([dir2save filesep 'resSummary_movie_%03d.mat'],iM));
%         for iM = 1 : numMovies
%             file2savePerMovie = fullfile(dir2save,sprintf('resSummary_movie_%03d.mat',iM));
%             if(redoPerMovieAnalysis || ~exist(file2savePerMovie,'file'))
%                 
% %                 resSummary(iM) = resultsIndTimeCoursePerMovie(ML.movieDataFile_{iM},file2savePerMovie);
% %             else
% %                 saved = load(file2savePerMovie);
% %                 resSummary(iM) = saved.resSummary;
%             end
%         end
        caseTimeList = [];
        caseResSummary = [];
        return;
    case 'load_only'
        for iM = 1 : numMovies
            file2savePerMovie = fullfile(dir2save,sprintf([dir2save filesep 'resSummary_movie_%03d.mat'],iM));
            % Block until file exists
            while(~exist(file2savePerMovie,'file'))
                % Pause for 5 seconds if analysis file does not exist
                disp(['Waiting 5 seconds for movie ' num2str(iM) ' out of ' num2str(numMovies)]);
                pause(5);
            end
            saved = load(file2savePerMovie);
            resSummary(iM) = saved.resSummary;
        end
    otherwise
        error(['Invalid parallel parameter: ' parallel]);
end

%% Sorting and Saving


if(~isempty(caseParam))
    %go over each case and put its results together
    for iCase = 1 : numCases

        %get parameters from input
        caseTimeList = caseParam(iCase).timeList;
        caseIndx = caseParam(iCase).indx;
        caseName = caseParam(iCase).name;
        caseMin0 = caseParam(iCase).indx0min;

        %sort time and add column for relative time
        offset = caseTimeList(caseMin0);
        [caseTimeList,indxSort] = sort(caseTimeList);
        caseTimeList = [caseTimeList caseTimeList-offset]; %#ok<AGROW>

        %collect and sort results
        caseResSummary = resSummary(caseIndx);
        caseResSummary = caseResSummary(indxSort);

        %Saves by default
        if saveFile

            %name variables properly for saving
            eval(['timeList_' caseName ' = caseTimeList;'])
            eval(['resSummary_' caseName ' = caseResSummary;'])

            %save results
            file2save = fullfile(dir2save,['resSummary_' caseName]); %#ok<NASGU>
            eval(['save(file2save,''timeList_' caseName ''',''resSummary_' caseName ''');']);

        end

    end
else
    caseTimeList = [];
    caseResSummary = resSummary;
    caseName = 'default';

    if saveFile
    
        %name variables properly for saving
        eval(['timeList_' caseName ' = caseTimeList;'])
        eval(['resSummary_' caseName ' = caseResSummary;'])
        
        %save results
        file2save = fullfile(dir2save,['resSummary_' caseName]); %#ok<NASGU>
        eval(['save(file2save,''timeList_' caseName ''',''resSummary_' caseName ''');']);
        
    end
end

