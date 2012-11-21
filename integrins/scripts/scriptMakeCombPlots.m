condName = {'AlphaV717Trunc'};
numMode = 4;
minNP = 4;

for iCond = 1 : length(condName);
    
    eval(['sptPropInWindow = sptPropInWindow' condName{iCond} 'CombModRatio;'])
    eval(['windowDistFromEdge = windowDistFromEdge' condName{iCond} 'Comb;'])
    
    for iType = 1
        
        %particle density
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.densityParticles,...
            windowDistFromEdge(iType),1,minNP(iCond),[condName{iCond} ', Particle density overall, protType ' num2str(iType)],...
            [condName{iCond} ' DensityOverall_protType' num2str(iType) '.fig'],[0.111 10 1],...
            'Density overall (ratio to onset)'); %,[-30 100 0 2.5]);
        
        %overall diffusion coefficient
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
            windowDistFromEdge(iType),1,minNP(iCond),[condName{iCond} ', Diffusion coefficient overall, protType ' num2str(iType)],...
            [condName{iCond} ' DiffCoefOverall_protType' num2str(iType) '.fig'],[0.111 10 1],...
            'Diffusion coefficient overall (ratio to onset)'); %,[-30 100 0.8 1.4]);
        
        %fraction in each diffusion mode
        for iMode = 1 : numMode(iCond)
            plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.fracModeClass,...
                windowDistFromEdge(iType),iMode,minNP(iCond),[condName{iCond} ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
                [condName{iCond} ' FractionMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1],...
                ['Fraction Mode ' num2str(iMode) ' (ratio to onset)']); %,[-30 100 0.1 1.4]);
        end
        
        %density of each diffusion mode
        for iMode = 1 : numMode(iCond)
            plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.densityModeClass,...
                windowDistFromEdge(iType),iMode,minNP(iCond),[condName{iCond} ', Density Mode ' num2str(iMode) ', protType ' num2str(iType)],...
                [condName{iCond} ' DensityMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1],...
                ['Density Mode ' num2str(iMode) ' (ratio to onset)']); %,[-30 100 0 10]);
        end
        
    end
    
end
