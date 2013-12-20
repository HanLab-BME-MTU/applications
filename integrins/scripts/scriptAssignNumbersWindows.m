
lengthMinMax = [5 99];

for i = 1 : length(movieStructAlphaVY773A)
    
    disp(num2str(i))
    
    activityLevel = movieStructAlphaVY773A(i).activityLevel;
    
    if activityLevel > 1
        
        tmp = movieStructAlphaVY773A(i).fileName{1};
        %         tmp = tmp(11:end);
        %         cd([tmp '/analysisAlphaVY773A/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis\adaptiveWindowsSym'])
        %         %%% for randomization test
        %         cd([topDir '\analysisAlphaVY773A\furtherAnalysis\randomizationTest\adaptiveWindows'])
        %         %%% for randomization test
        
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
        
        %         try
        %
        %             load windowNumbersAssignExt.mat
        %             windowNumbersAssignExtIn = windowNumbersAssignExt;
        %
        %             windowNumbersAssignExt = assignNumbers2Windows(tracksFinal,diffAnalysisRes,...
        %                 diffModeAnalysisRes,trackChar,windowsAll,1:400:maxFrame+1,...
        %                 windowTrackAssignExt,lengthMinMax,[2 3],windowNumbersAssignExtIn);
        %
        %         catch
        
        windowNumbersAssignExt = assignNumbers2Windows(tracksFinal,diffAnalysisRes,...
            diffModeAnalysisRes,trackChar,windowsAll,1:400:maxFrame+1,...
            windowTrackAssignExt,lengthMinMax);
        
        %         end
        
        Particles = windowNumbersAssignExt.Particles;
        save('windowNumbersAssignExt_Particles','Particles');
        NetDispPos = windowNumbersAssignExt.NetDispPos;
        save('windowNumbersAssignExt_NetDispPos','NetDispPos');
        NetDispNeg = windowNumbersAssignExt.NetDispNeg;
        save('windowNumbersAssignExt_NetDispNeg','NetDispNeg');
        Unclass = windowNumbersAssignExt.Unclass;
        save('windowNumbersAssignExt_Unclass','Unclass');
        Lin = windowNumbersAssignExt.Lin;
        save('windowNumbersAssignExt_Lin','Lin');
        Iso = windowNumbersAssignExt.Iso;
        save('windowNumbersAssignExt_Iso','Iso');
        IsoUnclass = windowNumbersAssignExt.IsoUnclass;
        save('windowNumbersAssignExt_IsoUnclass','IsoUnclass');
        Conf = windowNumbersAssignExt.Conf;
        save('windowNumbersAssignExt_Conf','Conf');
        Brown = windowNumbersAssignExt.Brown;
        save('windowNumbersAssignExt_Brown','Brown');
        Dir = windowNumbersAssignExt.Dir;
        save('windowNumbersAssignExt_Dir','Dir');
        ModeClass = windowNumbersAssignExt.ModeClass;
        save('windowNumbersAssignExt_ModeClass','ModeClass');
        Merge = windowNumbersAssignExt.Merge;
        save('windowNumbersAssignExt_Merge','Merge');
        Split = windowNumbersAssignExt.Split;
        save('windowNumbersAssignExt_Split','Split');
        
    end
    
end
