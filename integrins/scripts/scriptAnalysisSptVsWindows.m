
lengthMinMax = [5 99];

for i = 5 : length(movieStructAlphaV717Trunc)
    
    disp(num2str(i))
    
    activityLevel = movieStructAlphaV717Trunc(i).activityLevel;
    
    if activityLevel > 1
        
        tmp = movieStructAlphaV717Trunc(i).fileName{1};
        tmp = tmp(11:end);
        cd([tmp '/analysisAlphaV717Trunc/furtherAnalysis/adaptiveWindows'])
        
        sliceRange = movieStructAlphaV717Trunc(i).sliceRange;
        frameRange = movieStructAlphaV717Trunc(i).frameRange;
        
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
        
        maskDir = [tmp '/analysisCellEdgeModSmall/refined_masks/refined_masks_for_channel_1/'];
        tmp = dir([maskDir '*.tif']);
        firstMaskFile = [maskDir tmp(1).name];
        
        [sptPropInWindow,windowDistFromEdge,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
            tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
            protSamples,windowTrackAssignExt,windowNumbersAssignExt,...
            lengthMinMax,sliceRange,frameRange,firstMaskFile);
        
        save('particleBehaviorAdaptiveWindows120907','sptPropInWindow',...
            'windowDistFromEdge','analysisParam');
        
    end
    
end
