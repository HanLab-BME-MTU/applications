
for i = 1 : length(movieStructBeta3)

    disp(num2str(i))
    
    activityLevel = movieStructBeta3(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructBeta3(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisBeta3/furtherAnalysis'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisBeta3\furtherAnalysis'])
        %         %%% for randomization test
        %         cd([topDir '\analysisBeta3\furtherAnalysis\randomizationTest'])
        %         %%% for randomization test
                
        load tracksDiffusionLength5InMask.mat
        %         %%% for randomization test
        %         load tracksLength5InMaskRandom.mat
        %         %%% for randomization test
        
        cd adaptiveWindowsSym
        try
            load windowsActivityTracks.mat
        catch
            load windowsActivityTracks1.mat
            load windowsActivityTracks2.mat
        end
        
        trackChar = trackMotionCharProtrusion(tracksFinal,protSamples,trackWindowAssignComp,5);
        save('directTrackChar','trackChar');
        
    end
    
end

