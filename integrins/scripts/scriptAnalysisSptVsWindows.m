
lengthMinMax = [5 99];

for i = 1 : length(movieStructLifeact)
    
    disp(num2str(i))
    
    activityLevel = movieStructLifeact(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructLifeact(i).fileName{1};
        topDir = topDir(11:end);
        cd([topDir '/analysisLifeact/furtherAnalysis/adaptiveWindows'])
        
        sliceRange = movieStructLifeact(i).sliceRange;
        frameRange = movieStructLifeact(i).frameRange;
        
        load ../tracksDiffusionLength5InMask.mat
        load ../diffusionModeClassification.mat
        
        seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
        maxFrame = max(seqOfEvents(:,1));
        
        try
            load windowsActivityTracks.mat
        catch
            load windowsActivityTracks1.mat
            load windowsActivityTracks2.mat
            load windowsActivityTracks3.mat
            load windowsActivityTracks4.mat
            windowTrackAssignExt = cat(4,windowTrackAssignExt1,windowTrackAssignExt2);
            clear windowTrackAssignExt1 windowTrackAssignExt2
        end
        
        load directTrackChar.mat
        
        load windowNumbersAssignExt.mat
        
        maskDir = [topDir '/analysisCellEdgeSmall/SegmentationPackage/refined_masks/refined_masks_for_channel_1/'];
        tmp = dir([maskDir '*.tif']);
        firstMaskFile = [maskDir tmp(1).name];
        
        protWinParam = struct('numTypeProt',9,...
            'numPixInBand',2,...
            'numSmallBandsInBigBand',3,...
            'numBigBands',10,...
            'maxNegInc',10,...
            'maxPosInc',length(tmp)-1);
        
        edgePosStd = 1;
        
        [sptPropInWindow,windowDistFromEdge,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
            tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
            protSamples,windowTrackAssignExt,windowNumbersAssignExt,...
            lengthMinMax,sliceRange,frameRange,firstMaskFile,protWinParam,edgePosStd);
        
        save('particleBehaviorAdaptiveWindows121112','sptPropInWindow',...
            'windowDistFromEdge','analysisParam');
        
    end
    
end
