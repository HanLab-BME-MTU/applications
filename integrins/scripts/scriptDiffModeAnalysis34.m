
for i = 1 : length(movieStructLifeactFixed);
    
    disp(num2str(i));
    
    activityLevel = movieStructLifeactFixed(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructLifeactFixed(i).fileName;
        cd([tmp{1} '/analysisLifeact/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        
        %     [modeParam0,expParam0] = getDiffModes(tracksFinal,5,0.01,0,10,1,'test');
        %     [modeParam3,expParam3] = getDiffModes(tracksFinal,5,1,0,3,1,'test');
        %     [modeParam4,expParam4] = getDiffModes(tracksFinal,5,1,0,4,1,'test');
        %
        %     save('diffusionModeAnalysis34_1','modeParam0','expParam0','modeParam3',...
        %         'expParam3','modeParam4','expParam4');
        %
        %     movieStructLifeactFixed(i).diffModeParam0_1 = modeParam0;
        %     movieStructLifeactFixed(i).diffModeNum_1 = size(modeParam0,1);
        %     movieStructLifeactFixed(i).diffModeParam3_1 = modeParam3;
        %     movieStructLifeactFixed(i).diffModeParam4_1 = modeParam4;
        
        [modeParam0,numMode0,modeParamControl,numModeControl] = getDiffModes(tracksFinal,5,0.01,0,10,2,'test',[],1);
        %         modeParam3 = getDiffModes(tracksFinal,5,1,0,3,2,'test',[],0);
        %         modeParam4 = getDiffModes(tracksFinal,5,1,0,4,2,'test',[],0);
        modeParam3 = [];
        modeParam4 = [];
        
        save('diffusionModeAnalysis34_2New','modeParam0','numMode0','modeParam3',...
            'modeParam4','modeParamControl','numModeControl');
        
        diffModeAnalysisLifeactFixed(i).paramFree = modeParam0;
        diffModeAnalysisLifeactFixed(i).numMode = numMode0;
        diffModeAnalysisLifeactFixed(i).paramForced3 = modeParam3;
        diffModeAnalysisLifeactFixed(i).paramForced4 = modeParam4;
        diffModeAnalysisLifeactFixed(i).paramControlFree = modeParamControl;
        diffModeAnalysisLifeactFixed(i).numModeControl = numModeControl;
        
    end
    
end

