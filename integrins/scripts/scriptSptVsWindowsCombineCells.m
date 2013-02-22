
j = 0;

for i = 1 : length(movieStructFarnSim20p5)
    
    activityLevel = movieStructFarnSim20p5(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructFarnSim20p5(i).fileName{1};
        %         tmp = tmp(11:end);
        %         cd([tmp '/analysisFarnSim20p5/furtherAnalysis/adaptiveWindows'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisFarnSim20p5\furtherAnalysis\adaptiveWindowsSym'])
        
        %         %%% for randomization test
        %         cd([topDir '\analysisFarnSim20p5\furtherAnalysis\randomizationTest\adaptiveWindows'])
        %         %%% for randomization test
                
        load particleBehaviorAdaptiveWindows130219.mat
        sptPropInWindowFarnSim20p5Ind(:,j) = sptPropInWindow;
        windowDistFromEdgeFarnSim20p5Ind(:,j) = windowDistFromEdge;
        
    end
    
end

cd C:\kjData\Galbraiths\analysis\130219_sptVsWindowsCombCellsSym