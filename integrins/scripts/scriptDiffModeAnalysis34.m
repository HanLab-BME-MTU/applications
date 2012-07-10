
for i = 1 : length(analysisAlphaVInd);
    
    tmp = analysisAlphaVInd(i).fileName;
    cd([tmp{1} '/analysisAlphaV/furtherAnalysis'])
    
    load tracksDiffusionLength5InMask.mat
    
    [modeParam0,expParam0] = getDiffModes(tracksFinal,5,0.01,0,10,1,'test');
    [modeParam3,expParam3] = getDiffModes(tracksFinal,5,1,0,3,1,'test');
    [modeParam4,expParam4] = getDiffModes(tracksFinal,5,1,0,4,1,'test');
    
    save('diffusionModeAnalysis34_1','modeParam0','expParam0','modeParam3',...
        'expParam3','modeParam4','expParam4');
    
    analysisAlphaVInd(i).diffModeParam0_1 = modeParam0;
    analysisAlphaVInd(i).diffModeNum_1 = size(modeParam0,1);
    analysisAlphaVInd(i).diffModeParam3_1 = modeParam3;
    analysisAlphaVInd(i).diffModeParam4_1 = modeParam4;
    
    [modeParam0,expParam0] = getDiffModes(tracksFinal,5,0.01,0,10,2,'test');
    [modeParam3,expParam3] = getDiffModes(tracksFinal,5,1,0,3,2,'test');
    [modeParam4,expParam4] = getDiffModes(tracksFinal,5,1,0,4,2,'test');
    
    save('diffusionModeAnalysis34_2','modeParam0','expParam0','modeParam3',...
        'expParam3','modeParam4','expParam4');
    
    analysisAlphaVInd(i).diffModeParam0_2 = modeParam0;
    analysisAlphaVInd(i).diffModeNum_2 = size(modeParam0,1);
    analysisAlphaVInd(i).diffModeParam3_2 = modeParam3;
    analysisAlphaVInd(i).diffModeParam4_2 = modeParam4;
    
end
