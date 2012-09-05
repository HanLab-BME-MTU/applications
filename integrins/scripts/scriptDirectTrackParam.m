
for i = 1 : length(analysisLifeactInd);
    
    disp(num2str(i))
    
    tmp = analysisLifeactInd(i).fileName;
    cd([tmp{1} '/analysisLifeact/furtherAnalysis'])
    
    delete directTrackChar.mat
    
    activityLevel = analysisLifeactInd(i).activityLevel;
    
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

