
for i = 1 : length(movieStructFarnSim20p5)

    disp(num2str(i))
    
    activityLevel = movieStructFarnSim20p5(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructFarnSim20p5(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisFarnSim20p5/furtherAnalysis'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisFarnSim20p5\furtherAnalysis'])
        
        %         %%% for randomization test
        %         cd([topDir '\analysisFarnSim20p5\furtherAnalysis\randomizationTest'])
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

