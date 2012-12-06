
for i = 1 : length(movieStructAlphaVFixed);
    
    disp(num2str(i));
    
    activityLevel = movieStructAlphaVFixed(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructAlphaVFixed(i).fileName;
        cd([tmp{1} '/analysisAlphaV/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        
        [modeParam0,numMode0,modeParamControl,numModeControl] = getDiffModes(tracksFinal,5,0.01,0,10,2,'test',[],1);
        %         modeParam3 = getDiffModes(tracksFinal,5,1,0,3,2,'test',[],0);
        %         modeParam4 = getDiffModes(tracksFinal,5,1,0,4,2,'test',[],0);
        modeParam3 = [];
        modeParam4 = [];
        
        save('diffusionModeAnalysis34_2New','modeParam0','numMode0','modeParam3',...
            'modeParam4','modeParamControl','numModeControl');
        
        %         load diffusionModeAnalysis34_2New
        
        diffModeAnalysisAlphaVFixed(i).paramFree = modeParam0;
        diffModeAnalysisAlphaVFixed(i).numMode = numMode0;
        diffModeAnalysisAlphaVFixed(i).paramForced3 = modeParam3;
        diffModeAnalysisAlphaVFixed(i).paramForced4 = modeParam4;
        diffModeAnalysisAlphaVFixed(i).paramControlFree = modeParamControl;
        diffModeAnalysisAlphaVFixed(i).numModeControl = numModeControl;
        
    end
    
end

