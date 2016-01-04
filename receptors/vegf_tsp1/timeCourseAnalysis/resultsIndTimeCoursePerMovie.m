function [ resSummary ] = resultsIndTimeCoursePerMovie( MD, saveFile, channels )
%resultsIndTimeCoursePerMovie Subroutine of resultsIndTimeCourse
%
% INPUT
% MD        : MovieData object or string for MovieData .mat file
% saveFile  : File to which to save resSummary (optional, default: no save)
% channels : (optional) What channels to use, default:
%                       1:length(MD.channels_)
% See also resultsIndTimeCourse

% Based on resultsIndTimeCourse by Khuloud Jaqaman, March 2015
% Mark Kittisopikul, October2015

if(ischar(MD))
    MD = MovieData.load(MD);
end
if(nargin < 2 || isempty(channels))   
    channels = 1 : length(MD.channels_);
    curChannels = channels;
else
    curChannels = channels(:)';
    curChannels = curChannels(curChannels <= length(MD.channels_));
end
   
    % channels could be empty or exceed the number of channels in this movie
    resSummary(1,max(channels)) = timeCourseAnalysis.util.emptyResSummaryStruct;

    progressTextMultiple('channel', length(curChannels));

    
    for iC = curChannels

        iProcDiff = MD.getProcessIndex('MotionAnalysisProcess',1,0); %diffusion analysis and tracks
        if isempty(iProcDiff)
            error([MD.movieDataPath_ ' : Process Missing']);
        end
        motionAnalysis = load(MD.processes_{iProcDiff}.outFilePaths_{iC});
        iProcMask = MD.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
        if ~isempty(iProcMask)
            % TODO: Is there a different mask for the channel?
            mask = imread(fullfile(MD.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
            cellArea = sum(mask(:));
        else
            mask = [];
            cellArea = 512^2; %hard-code for now - look for better solution given the peculiarities of our simultaneous 2-color imaging
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
        %KJ, 151121: get amplitudes only in first 20 frames of movie
        %(instead of all throughout), to minimize effect of photobleaching.
        %This means not all tracks will get an amplitude
        ampMeanPerTraj = nanmean(ampMat(:,1:20),2);

        %average properties per motion class
        [ampMeanPerClass,diffCoefMeanPerClass,confRadMeanPerClass] = deal(NaN(5,1));
        for i = 1 : 2
            indxClass = find(trajClass == (i-1));
            ampMeanPerClass(i) = nanmean(ampMeanPerTraj(indxClass));
            diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
            confRadMeanPerClass(i) = mean(trajConfRad(indxClass));
        end
        for i = 3 : 4
            indxClass = find(trajClass == (i-1));
            ampMeanPerClass(i) = nanmean(ampMeanPerTraj(indxClass));
            diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
        end
        ampMeanPerClass(5) = nanmean(ampMeanPerTraj(isnan(trajClass)));

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

        %normalize mean amplitudes by first mode mean
        %for amp per class and first 20 frames, use mode of first 20 frames
        ampMeanPerClass = [ampMeanPerClass ampMeanPerClass/ampStatsF20(2)]; %#ok<AGROW>
        ampStatsF20 = [ampStatsF20 ampStatsF20(1)/ampStatsF20(2)]; %#ok<AGROW>
        %for last 20 frames, use mode of last 20 frames
        ampStatsL20 = [ampStatsL20 ampStatsL20(1)/ampStatsL20(2)]; %#ok<AGROW>

        %merge and split statistics
        statsMS = calcStatsMS_noMotionInfo(motionAnalysis.tracks,5,1,1);
        tmp = calcMergeSplitTimes(motionAnalysis.tracks,5,[],1);
        tmp = tmp.all;
        msTimeInfo = [mean(tmp.timeMerge2Split) mean(tmp.timeSplit2MergeSelf) mean(tmp.timeSplit2MergeOther) mean(tmp.timeMerge2End) mean(tmp.timeStart2Split)];
        %         msTimeInfo = tmp.timeMerge2Split;

        %results for output
        resSummary(1,iC).diffSummary = diffSummary;
        resSummary(1,iC).diffCoefMeanPerClass = diffCoefMeanPerClass;
        resSummary(1,iC).confRadMeanPerClass = confRadMeanPerClass;
        resSummary(1,iC).ampMeanPerClass = ampMeanPerClass;
        resSummary(1,iC).ampStatsF20 = ampStatsF20;
        resSummary(1,iC).ampStatsL20 = ampStatsL20;
        resSummary(1,iC).statsMS = statsMS;
        resSummary(1,iC).msTimeInfo = msTimeInfo;
        resSummary(1,iC).cellArea = cellArea;
        
        %Progress Counter
        progressTextMultiple();
    end
    
    if(nargin > 1 && saveFile)
        try
            save(saveFile,'resSummary');
        catch err
            warning(['resultIndTimeCoursePerMovie: Could not save to file, ' saveFile]);
            disp(getReport(err));
        end
    end
    
    progressTextMultiple();

end

