
for i = 1 : length(movieStructAlphaVPax);
    
    disp(num2str(i))
    
    if movieStructAlphaVPax(i).activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName{1};
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaV\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        spatialRadius = 4;
        timeRadius = 50;
        numNeighbors = numNeighborsTrack(tracksFinal,spatialRadius,timeRadius);
        save('numNeighborsPerTrack','numNeighbors','spatialRadius','timeRadius');
        
    end
    
end
