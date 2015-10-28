function [caseResSummary] = resultsIndTimeCourseMod(ML, saveFile)
%RESULTSINDTIMECOURSE compiles the results of a group of movies making one or multiple timecourse datasets and orders them based on time
%
%SYNOPSIS [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam)
%
%INPUT  ML       : MovieList object containing all movies, either all
%                  belonging to one timecourse.
%        saveFile: Logical that determines if this function saves a file.
%                  The default is 'true'. 0 or 1 instead of true or false
%                  will work.
%    
%OUTPUT 
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
%modified from resultsIndTimeCourse
%Tae Kim, July 2015

%% Input and pre-processing





if nargin < 2
    saveFile = true;
end

%get number of movies and number of cases
numMovies = length(ML.movieDataFile_);
numCases = 1;

%reserve memory for individual movie results
resSummary = struct('diffSummary',[],'diffCoefMeanPerClass',[],...
    'confRadMeanPerClass',[],'ampMeanPerClass',[],'ampStatsF20',[],...
    'ampStatsL20',[],'statsMS',[]);

resSummary = repmat(resSummary,numMovies,1);

%define directory for saving results
dir2save = [ML.movieListPath_ filesep 'analysisKJ'];

%% Calculations
nMD = numMovies;
iMD = 0;
%printLength = fprintf(1,'%g/%g MovieData analyzed\n', iMD, nMD);
progressTextMultiple('analyzing MD', nMD)

%go over all movies
for iM = 1 : numMovies
    
    %read in raw results
    %     disp([num2str(iM) '  ' ML.movieDataFile_{iM}]);
    %No need to load MD again if it has already been loaded
    if isempty(ML.movies_{iM})
        MD = MovieData.load(ML.movieDataFile_{iM});
    else
        MD = ML.movies_{iM};
    end
    
    iProcDiff = MD.getProcessIndex('MotionAnalysisProcess',1,0); %diffusion analysis and tracks
    if isempty(iProcDiff)
        error([MD.movieDataPath_ ' : Process Missing']);
    end
    load(MD.processes_{iProcDiff}.outFilePaths_{1});
    iProcMask = MD.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
    if ~isempty(iProcMask)
        mask = imread(fullfile(MD.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
    else
        mask = [];
    end
    
    %limit analysis to tracks in mask if supplied
    if ~isempty(mask) && any(mask(:)==0)
        %keep only tracks in mask
        numTracks = length(tracks);
        keepTrack = ones(numTracks,1);
        for iTrack = 1 : numTracks
            xCoord = tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
            yCoord = tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
            meanPosX = round(nanmean(xCoord(:)));
            meanPosY = round(nanmean(yCoord(:)));
            keepTrack(iTrack) = mask(meanPosY,meanPosX);
        end
        indxKeep = find(keepTrack);
        tracks = tracks(indxKeep);
        diffAnalysisRes = diffAnalysisRes(indxKeep);
        %redo diffusion analysis summary
        if length(indxKeep) < numTracks
            minTrackLen = 5;
            probDim = 2;
            extractType = 1;
            [probMotionType,motionChar] = summarizeDiffAnRes(tracks,minTrackLen,probDim,diffAnalysisRes,extractType);
            diffAnalysisRes(1).summary.probMotionType = probMotionType;
            diffAnalysisRes(1).summary.motionChar = motionChar;
        end
    end
    
    %diffusion analysis summary
    diffSummary = diffAnalysisRes(1).summary;
    
    %diffusion analysis per track
    trajClass = vertcat(diffAnalysisRes.classification);
    trajClass = trajClass(:,2);
    trajDiffCoef = catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef');
    trajConfRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');
    
    %amplitude matrix
    tracksMat = convStruct2MatIgnoreMS(tracks);
    ampMat = tracksMat(:,4:8:end);
    
    %amplitude per track
    ampMeanPerTraj = nanmean(ampMat,2);
    
    %average properties per motion class
    [ampMeanPerClass,diffCoefMeanPerClass,confRadMeanPerClass] = deal(NaN(5,1));
    for i = 1 : 2
        indxClass = find(trajClass == (i-1));
        ampMeanPerClass(i) = mean(ampMeanPerTraj(indxClass));
        diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
        confRadMeanPerClass(i) = mean(trajConfRad(indxClass));
    end
    for i = 3 : 4
        indxClass = find(trajClass == (i-1));
        ampMeanPerClass(i) = mean(ampMeanPerTraj(indxClass));
        diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
    end
    ampMeanPerClass(5) = mean(ampMeanPerTraj(isnan(trajClass)));
    
    %amplitude statistics in first 20 frames
    ampVec = ampMat(:,1:20);
    ampVec = ampVec(~isnan(ampVec));
    %tic;
    [~,~,modeParam] = fitHistWithGaussians(ampVec,0.01,0,3,0,[],2,[],1,[],0);
    %toc;%--------------------------------------------------------------------------------------------------------------------------------------
    numMode = size(modeParam,1);
    modeParamMean = modeParam(1,1);
    modeParamStd  = modeParam(1,2);
    ampMode1Mean = exp( modeParamMean + modeParamStd.^2/2 );
    ampMode1Std  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
    ampMode1Frac = modeParam(1,4)/sum(modeParam(:,4));
    ampStatsF20 = [mean(ampVec) ampMode1Mean ampMode1Std ampMode1Frac numMode];
            
    %amplitude statistics in last 20 frames
    ampVec = ampMat(:,end-19:end);
    ampVec = ampVec(~isnan(ampVec));
    %tic;
    [~,~,modeParam] = fitHistWithGaussians(ampVec,0.01,0,3,0,[],2,[],1,[],0); %toc;%-----------------------------------------------------
    numMode = size(modeParam,1);
    modeParamMean = modeParam(1,1);
    modeParamStd  = modeParam(1,2);
    ampMode1Mean = exp( modeParamMean + modeParamStd.^2/2 );
    ampMode1Std  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
    ampMode1Frac = modeParam(1,4)/sum(modeParam(:,4));
    ampStatsL20 = [mean(ampVec) ampMode1Mean ampMode1Std ampMode1Frac numMode];
    
    %normalize mean amplitudes by first mode mean from last 20 frames
    ampMeanPerClass = [ampMeanPerClass ampMeanPerClass/ampStatsL20(2)]; %#ok<AGROW>
    ampStatsF20 = [ampStatsF20 ampStatsF20(1)/ampStatsL20(2)]; %#ok<AGROW>
    ampStatsL20 = [ampStatsL20 ampStatsL20(1)/ampStatsL20(2)]; %#ok<AGROW>

    %merge and split statistics
    statsMS = calcStatsMS_noMotionInfo(tracks,5,1,1);
    
    %results for output
    resSummary(iM).diffSummary = diffSummary;
    resSummary(iM).diffCoefMeanPerClass = diffCoefMeanPerClass;
    resSummary(iM).confRadMeanPerClass = confRadMeanPerClass;
    resSummary(iM).ampMeanPerClass = ampMeanPerClass;
    resSummary(iM).ampStatsF20 = ampStatsF20;
    resSummary(iM).ampStatsL20 = ampStatsL20;
    resSummary(iM).statsMS = statsMS;
    
    %Progress Counter
    progressTextMultiple();
    %{
    fprintf(repmat('\b',1,printLength));
    iMD = iMD + 1;
    printLength = fprintf(1,'%g/%g MovieData analyzed\n', iMD, nMD);
    %}
end
%fprintf(repmat('\b',1,printLength));

%% Sorting and Saving

%go over each case and put its results together
for iCase = 1 : numCases
    

    
    %collect and sort results
    caseResSummary = resSummary;
    
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

