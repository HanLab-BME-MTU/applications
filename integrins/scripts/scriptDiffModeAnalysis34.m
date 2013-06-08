
for i = 1 : length(movieStructBeta3);
    
    disp(num2str(i));
    
    activityLevel = movieStructBeta3(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructBeta3(i).fileName{1};
        %         cd([tmp '/analysisAlphaV/furtherAnalysis'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisBeta3\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        
        [modeParam0,numMode0,modeParamControl,numModeControl] = getDiffModes(tracksFinal,5,0.01,0,10,2,'test',[],1);
        %         modeParam3 = getDiffModes(tracksFinal,5,1,0,3,2,'test',[],0);
        %         modeParam4 = getDiffModes(tracksFinal,5,1,0,4,2,'test',[],0);
        modeParam3 = [];
        modeParam4 = [];
        
        save('diffusionModeAnalysis34_2New','modeParam0','numMode0','modeParam3',...
            'modeParam4','modeParamControl','numModeControl');
        
        %         load diffusionModeAnalysis34_2New
        
        diffModeAnalysisBeta3(i).paramFree = modeParam0;
        diffModeAnalysisBeta3(i).numMode = numMode0;
        diffModeAnalysisBeta3(i).paramForced3 = modeParam3;
        diffModeAnalysisBeta3(i).paramForced4 = modeParam4;
        diffModeAnalysisBeta3(i).paramControlFree = modeParamControl;
        diffModeAnalysisBeta3(i).numModeControl = numModeControl;
        
    end
    
end

