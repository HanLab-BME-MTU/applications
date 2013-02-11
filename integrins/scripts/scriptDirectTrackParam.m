
for i = 1 : length(movieStructFarn)

    disp(num2str(i))
    
    activityLevel = movieStructFarn(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructFarn(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisFarn/furtherAnalysis'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        %         cd([topDir '\analysisFarn\furtherAnalysis'])
        
        %%% for randomization test
        cd([topDir '\analysisFarn\furtherAnalysis\randomizationTest'])
        %%% for randomization test
                
        %         load tracksDiffusionLength5InMask.mat
        
        %%% for randomization test
        load tracksLength5InMaskRandom.mat
        %%% for randomization test
        
        cd adaptiveWindows
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

