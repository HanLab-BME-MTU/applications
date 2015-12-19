function [ resSummary ] = resultsIndTimeCoursePerMovie( MD, saveFile )
%resultsIndTimeCoursePerMovie Subroutine of resultsIndTimeCourse
%
% INPUT
% MD        : MovieData object or string for MovieData .mat file
% saveFile  : File to which to save resSummary (optional, default: no save)
%
% See also resultsIndTimeCourse

% Based on resultsIndTimeCourse by Khuloud Jaqaman, March 2015
% Mark Kittisopikul, October2015

if(ischar(MD))
    MD = MovieData.load(MD);
end

    iProcDiff = MD.getProcessIndex('MotionAnalysisProcess',1,0); %diffusion analysis and tracks
    motionAnalysis = load(MD.processes_{iProcDiff}.outFilePaths_{1});
    iProcMask = MD.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
    if ~isempty(iProcMask)
        mask = imread(fullfile(MD.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
    else
        mask = [];
    end
    
    %limit analysis to tracks in mask if supplied
    if ~isempty(mask) && any(mask(:)==0)
        %keep only tracks in mask
        numTracks = length(motionAnalysis.tracks);
        keepTrack = ones(numTracks,1);
        for iTrack = 1 : numTracks
            xCoord = motionAnalysis.tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
            yCoord = motionAnalysis.tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
            meanPosX = round(nanmean(xCoord(:)));
            meanPosY = round(nanmean(yCoord(:)));
            keepTrack(iTrack) = mask(meanPosY,meanPosX);
        end
        indxKeep = find(keepTrack);
        motionAnalysis.tracks = motionAnalysis.tracks(indxKeep);
        motionAnalysis.diffAnalysisRes = motionAnalysis.diffAnalysisRes(indxKeep);
        %redo diffusion analysis summary
        if length(indxKeep) < numTracks
            minTrackLen = 5;
            probDim = 2;
            extractType = 1;
            [probMotionType,motionChar] = summarizeDiffAnRes(motionAnalysis.tracks,minTrackLen,probDim,motionAnalysis.diffAnalysisRes,extractType);
            motionAnalysis.diffAnalysisRes(1).summary.probMotionType = probMotionType;
            motionAnalysis.diffAnalysisRes(1).summary.motionChar = motionChar;
        end
    end
    
    %diffusion analysis summary
    diffSummary = motionAnalysis.diffAnalysisRes(1).summary;
    
    %diffusion analysis per track
    trajClass = vertcat(motionAnalysis.diffAnalysisRes.classification);
    trajClass = trajClass(:,2);
    trajDiffCoef = catStruct(1,'motionAnalysis.diffAnalysisRes.fullDim.normDiffCoef');
    trajConfRad = catStruct(1,'motionAnalysis.diffAnalysisRes.confRadInfo.confRadius');
    
    %amplitude matrix
    tracksMat = convStruct2MatIgnoreMS(motionAnalysis.tracks);
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
    [~,~,modeParam] = fitHistWithGaussians(ampVec,0.05,0,3,0,[],2,[],1,[],0);
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
    [~,~,modeParam] = fitHistWithGaussians(ampVec,0.05,0,3,0,[],2,[],1,[],0);
    numMode = size(modeParam,1);
    modeParamMean = modeParam(1,1);
    modeParamStd  = modeParam(1,2);
    ampMode1Mean = exp( modeParamMean + modeParamStd.^2/2 );
    ampMode1Std  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
    ampMode1Frac = modeParam(1,4)/sum(modeParam(:,4));
    ampStatsL20 = [mean(ampVec) ampMode1Mean ampMode1Std ampMode1Frac numMode];
    
    %normalize mean amplitudes by first mode mean from last 20 frames
    ampMeanPerClass = [ampMeanPerClass ampMeanPerClass/ampStatsL20(2)]; 
    ampStatsF20 = [ampStatsF20 ampStatsF20(1)/ampStatsL20(2)]; 
    ampStatsL20 = [ampStatsL20 ampStatsL20(1)/ampStatsL20(2)]; 

    %merge and split statistics
    statsMS = calcStatsMS_noMotionInfo(motionAnalysis.tracks,5,1,1);
    
    %results for output
    resSummary.diffSummary = diffSummary;
    resSummary.diffCoefMeanPerClass = diffCoefMeanPerClass;
    resSummary.confRadMeanPerClass = confRadMeanPerClass;
    resSummary.ampMeanPerClass = ampMeanPerClass;
    resSummary.ampStatsF20 = ampStatsF20;
    resSummary.ampStatsL20 = ampStatsL20;
    resSummary.statsMS = statsMS;
    
    if(nargin > 1)
        try
            save(saveFile,'resSummary');
        catch err
            warning(['resultIndTimeCoursePerMovie: Could not save to file, ' saveFile]);
            disp(getReport(err));
        end
    end

end

