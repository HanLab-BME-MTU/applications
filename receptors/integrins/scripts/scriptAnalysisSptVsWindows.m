
lengthMinMax = [5 99];

for i = 12 : length(movieStructAlphaVY773A)
    
    disp(num2str(i))
    
    activityLevel = movieStructAlphaVY773A(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructAlphaVY773A(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisAlphaVY773A/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis\adaptiveWindowsSym'])
        %         %% for randomization test
        %         cd([topDir '\analysisAlphaVY773A\furtherAnalysis\randomizationTest\adaptiveWindowsSym'])
        %         %% for randomization test
        
        sliceRange = movieStructAlphaVY773A(i).sliceRange;
        frameRange = movieStructAlphaVY773A(i).frameRange;
        
        load ../tracksDiffusionLength5InMask.mat
        load ../diffusionModeClassification.mat
        %         %%% for randomization test
        %         load ../../tracksDiffusionLength5InMask.mat
        %         load ../tracksLength5InMaskRandom.mat
        %         load ../../diffusionModeClassification.mat
        %         %%% for randomization test
        
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
        
        try
            load windowNumbersAssignExt.mat
        catch
            load windowNumbersAssignExt_Particles.mat
            windowNumbersAssignExt.Particles = Particles;
            load windowNumbersAssignExt_NetDispPos.mat
            windowNumbersAssignExt.NetDispPos = NetDispPos;
            load windowNumbersAssignExt_NetDispNeg.mat
            windowNumbersAssignExt.NetDispNeg = NetDispNeg;
            load windowNumbersAssignExt_Unclass.mat
            windowNumbersAssignExt.Unclass = Unclass;
            load windowNumbersAssignExt_Lin.mat
            windowNumbersAssignExt.Lin = Lin;
            load windowNumbersAssignExt_Iso.mat
            windowNumbersAssignExt.Iso = Iso;
            load windowNumbersAssignExt_IsoUnclass.mat
            windowNumbersAssignExt.IsoUnclass = IsoUnclass;
            load windowNumbersAssignExt_Conf.mat
            windowNumbersAssignExt.Conf = Conf;
            load windowNumbersAssignExt_Brown.mat
            windowNumbersAssignExt.Brown = Brown;
            load windowNumbersAssignExt_Dir.mat
            windowNumbersAssignExt.Dir = Dir;
            load windowNumbersAssignExt_ModeClass.mat
            windowNumbersAssignExt.ModeClass = ModeClass;
            load windowNumbersAssignExt_Merge.mat
            windowNumbersAssignExt.Merge = Merge;
            load windowNumbersAssignExt_Split.mat
            windowNumbersAssignExt.Split = Split;
        end
        
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
        
        save('particleBehaviorAdaptiveWindows130326','sptPropInWindow',...
            'windowDistFromEdge','analysisParam');
        
    end
    
end
