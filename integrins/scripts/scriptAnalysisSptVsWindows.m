
lengthMinMax = [5 99];

% for i = 1 : length(movieStructAlphaV)
for i = 1 : 1
    
    disp(num2str(i))
    
    activityLevel = movieStructAlphaV(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructAlphaV(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisAlphaV/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        %         cd([topDir '\analysisAlphaV\furtherAnalysis\adaptiveWindows'])
        
        %%% for randomization test
        cd([topDir '\analysisAlphaV\furtherAnalysis\randomizationTest\adaptiveWindows'])
        %%% for randomization test
        
        sliceRange = movieStructAlphaV(i).sliceRange;
        frameRange = movieStructAlphaV(i).frameRange;
        
        %         load ../tracksDiffusionLength5InMask.mat
        %         load ../diffusionModeClassification.mat
        
        %%% for randomization test
        load ../../tracksDiffusionLength5InMask.mat
        load ../tracksLength5InMaskRandom.mat
        load ../../diffusionModeClassification.mat
        %%% for randomization test
        
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
        
        maxFrame = (size(windowsAll,1)-1)*400;

        load directTrackChar.mat
        
        load windowNumbersAssignExt.mat
        
        %         maskDir = [topDir '/analysisCellEdgeSmall/SegmentationPackage/refined_masks/refined_masks_for_channel_1/'];
        maskDir = [topDir '\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\'];
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
        
        save('particleBehaviorAdaptiveWindowsRandom130210','sptPropInWindow',...
            'windowDistFromEdge','analysisParam');
        
    end
    
end
