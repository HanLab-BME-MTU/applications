
j = 0;

for i = 1 : length(movieStructBeta3)
    
    activityLevel = movieStructBeta3(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructBeta3(i).fileName{1};
        %         tmp = tmp(11:end);
        %         cd([tmp '/analysisBeta3/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisBeta3\furtherAnalysis\adaptiveWindowsSym'])
        %         %%% for randomization test
        %         cd([topDir '\analysisBeta3\furtherAnalysis\randomizationTest\adaptiveWindowsSym'])
        %         %%% for randomization test
                
        load particleBehaviorAdaptiveWindows130326.mat
        sptPropInWindowBeta3Ind(:,j) = sptPropInWindow;
        windowDistFromEdgeBeta3Ind(:,j) = windowDistFromEdge;
        
    end
    
end

cd C:\kjData\Galbraiths\analysis\130327_sptVsWindowsSymCompArea
