
j = 0;

for i = 1 : length(movieStructAlphaVY773A)
    
    activityLevel = movieStructAlphaVY773A(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructAlphaVY773A(i).fileName{1};
        %         tmp = tmp(11:end);
        %         cd([tmp '/analysisAlphaVY773A/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis\adaptiveWindowsSym'])
        %         %%% for randomization test
        %         cd([topDir '\analysisAlphaVY773A\furtherAnalysis\randomizationTest\adaptiveWindowsSym'])
        %         %%% for randomization test
                
        load particleBehaviorAdaptiveWindows130326.mat
        sptPropInWindowAlphaVY773AInd(:,j) = sptPropInWindow;
        windowDistFromEdgeAlphaVY773AInd(:,j) = windowDistFromEdge;
        
    end
    
end

cd C:\kjData\Galbraiths\analysis\130327_sptVsWindowsSymCompArea
