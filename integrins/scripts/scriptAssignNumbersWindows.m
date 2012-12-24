
lengthMinMax = [5 99];

for i = 1 : length(movieStructLifeact)
    
    disp(num2str(i))
    
    activityLevel = movieStructLifeact(i).activityLevel;
    
    if activityLevel > 1 && i ~= 6
        
        tmp = movieStructLifeact(i).fileName{1};
        tmp = tmp(11:end);
        cd([tmp '/analysisLifeact/furtherAnalysis/adaptiveWindows'])
        
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
        
        save('windowNumbersAssignExt','windowNumbersAssignExt','-v7.3')
        
    end
    
end
