condName = {'Lifeact'}; %,'AlphaV717Trunc','AlphaV724Trunc','Beta3','Farn','Talin','Lifeact'};
% condName = {'AlphaVY773A','AlphaVY773ATmp'};
% condName = {'AlphaVSim20p5','AlphaVSim22','FarnSim20p5','FarnSim22'};
numMode  = [4 4 4 4 4 4 2];
minNP    = [5 5 5 5 5 5 5];

for iCond = 1 : length(condName);
    
    eval(['sptPropInWindow = sptPropInWindow' condName{iCond} 'CombModRatio;'])
    eval(['windowDistFromEdge = windowDistFromEdge' condName{iCond} 'Comb;'])
    
    for iType = 1
        
        %particle density
        testTrendSignificanceV3paper(sptPropInWindow(iType).directAll.densityParticles,...
            windowDistFromEdge(iType),1,minNP(iCond),[condName{iCond} ', Particle density overall, protType ' num2str(iType)],...
            [condName{iCond} ' DensityOverall_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],...
            'Molecule density (per um^2)',[-60 70 0 11],1,[condName{iCond} ', Edge & windows, protType ' num2str(iType)],...
            [condName{iCond} ' EdgeWindows_protType' num2str(iType) '.fig']);
        
        %overall diffusion coefficient
        testTrendSignificanceV3paper(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
            windowDistFromEdge(iType),1,minNP(iCond),[condName{iCond} ', Diffusion coefficient overall, protType ' num2str(iType)],...
            [condName{iCond} ' DiffCoefOverall_protType' num2str(iType) '.fig'],[0.111 10 (0.111^2)/0.025],...
            'Diffusion coefficient (um^2/s)',[-60 70 0.085 0.22]);
        
        %         %fraction in each diffusion mode
        %         for iMode = 1 : numMode(iCond)
        %             plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).diffModeAnalysis.fracModeClass,...
        %                 windowDistFromEdge(iType),iMode,minNP(iCond),[condName{iCond} ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
        %                 [condName{iCond} ' FractionMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1],...
        %                 ['Fraction Mode ' num2str(iMode)],[-100 150 0 1]);
        %         end
        %
        %         %density of each diffusion mode
        %         for iMode = 1 : numMode(iCond)
        %             plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).diffModeAnalysis.densityModeClass,...
        %                 windowDistFromEdge(iType),iMode,minNP(iCond),[condName{iCond} ', Density Modes '  num2str(iMode) ', protType ' num2str(iType)],...
        %                 [condName{iCond} ' DensityMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],...
        %                 ['Density Mode' num2str(iMode) ' (per um^2)'],[-100 150 0 6.5]);
        %         end
        
    end
    
end
