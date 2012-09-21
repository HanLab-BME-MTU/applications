
for i = 1 : length(movieStructAlphaV);
    
    disp(num2str(i))
    
    topDir = movieStructAlphaV(i).fileName{1};
    cd([topDir '/analysisAlphaV/furtherAnalysis'])
    
    %     delete directTrackChar.mat
    
    activityLevel = movieStructAlphaV(i).activityLevel;
    
    if activityLevel > 1
        
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

