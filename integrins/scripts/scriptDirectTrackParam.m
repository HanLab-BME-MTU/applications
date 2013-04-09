
for i = 1 : length(movieStructAlphaVY773A)

    disp(num2str(i))
    
    activityLevel = movieStructAlphaVY773A(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructAlphaVY773A(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisAlphaVY773A/furtherAnalysis'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis'])
        %         %%% for randomization test
        %         cd([topDir '\analysisAlphaVY773A\furtherAnalysis\randomizationTest'])
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

