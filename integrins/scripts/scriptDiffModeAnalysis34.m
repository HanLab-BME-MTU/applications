
for i = 1 : length(movieStructAlphaVY773A);
    
    disp(num2str(i));
    
    activityLevel = movieStructAlphaVY773A(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructAlphaVY773A(i).fileName{1};
        %         cd([tmp '/analysisAlphaV/furtherAnalysis'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        
        [modeParam0,numMode0,modeParamControl,numModeControl] = getDiffModes(tracksFinal,5,0.01,0,10,2,'test',[],1);
        %         modeParam3 = getDiffModes(tracksFinal,5,1,0,3,2,'test',[],0);
        %         modeParam4 = getDiffModes(tracksFinal,5,1,0,4,2,'test',[],0);
        modeParam3 = [];
        modeParam4 = [];
        
        save('diffusionModeAnalysis34_2New','modeParam0','numMode0','modeParam3',...
            'modeParam4','modeParamControl','numModeControl');
        
        %         load diffusionModeAnalysis34_2New
        
        diffModeAnalysisAlphaVY773A(i).paramFree = modeParam0;
        diffModeAnalysisAlphaVY773A(i).numMode = numMode0;
        diffModeAnalysisAlphaVY773A(i).paramForced3 = modeParam3;
        diffModeAnalysisAlphaVY773A(i).paramForced4 = modeParam4;
        diffModeAnalysisAlphaVY773A(i).paramControlFree = modeParamControl;
        diffModeAnalysisAlphaVY773A(i).numModeControl = numModeControl;
        
    end
    
end

