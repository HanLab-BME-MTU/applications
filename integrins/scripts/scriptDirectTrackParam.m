
for i = 1 : length(movieStructLifeact)
    
    disp(num2str(i))
    
    activityLevel = movieStructLifeact(i).activityLevel;
    
    if activityLevel > 1
        
        topDir = movieStructLifeact(i).fileName{1};
        topDir = topDir(11:end);
        cd([topDir '/analysisLifeact/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        
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

