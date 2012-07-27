clear all
close all

windowRange = (50:148)';
% windowRange = [];
frameRange = [];

windowsAll = putWindowsTogether;
load ../../../analysisCellEdgeModSmall/protrusion_samples/protrusion_samples.mat
load ../tracksDiffusionLength5InMask.mat

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[windowTrackAssign,trackWindowAssign,trackWindowAssignComp,windowTrackAssignExt] = ...
    assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);
save('windowsActivityTracks','protSamples','windowsAll','trackWindowAssign',...
    'trackWindowAssignComp','windowTrackAssign','windowTrackAssignExt')

load ../diffusionModeClassification.mat
load ../directTrackChar.mat

[sptPropInWindow,~,~,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
    tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
    1:400:maxFrame+1,protSamples,5,windowRange,frameRange,windowTrackAssignExt);

save('particleBehaviorAdaptiveWindows','sptPropInWindow','analysisParam');

for iType = [1 2 9]
    
    %particle density
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.spDensity,...
        20,'Overall particle density',['figuresSptVsWindows/spDensity_protType' num2str(iType) '.fig']);
    
    %frame-to-frame displacement
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.f2fDispMag2D,...
        20,'Frame-to-frame displacement',['figuresSptVsWindows/f2fDisp_protType' num2str(iType) '.fig']);
    
    %overall diffusion coefficient
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
        20,'Overall diffusion coefficient',['figuresSptVsWindows/diffCoef_protType' num2str(iType) '.fig']);
    
    %fraction in each diffusion mode
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.fracMode,...
        20,'Diffusion mode fraction',['figuresSptVsWindows/modeFraction_protType' num2str(iType) '.fig']);
    
    %particle density in each diffusion mode
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.modeDensity,...
        20,'Diffusion mode density',['figuresSptVsWindows/modeDensity_protType' num2str(iType) '.fig']);
    
    %diffusion coefficient in each diffusion mode
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeInd,...
        20,'Diffusion coefficient per mode',['figuresSptVsWindows/modeDiffCoef_protType' num2str(iType) '.fig']);
    
end


% windowRange = windowMotionChar.indxGoodWindow;
% frameRange = [];
%
% [sptPropInWindow,~,~,analysisParam] = ...
%     particleBehaviorRelToActivityOnset(tracksFinal,diffAnalysisRes,...
%     diffModeAnalysisRes,trackChar,windowsAll,1:400:maxFrame+1,protSamples,...
%     5,[1 1; 2 2; 3 3],windowRange,frameRange,windowTrackAssign);
% 
% save('particleBehaviorRelActivityGood','sptPropInWindow','analysisParam');
% 
% for iBand = 1 : 2
%     
%     %fraction in each diffusion mode
%     plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.fracMode,...
%         'Fraction in diffusion mode',['figuresSptVsWindows/fracModes_band' num2str(iBand) '_good.fig']);
%     
%     %overall diffusion coefficient
%     plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.diffCoefModeAll,...
%         'Overall diffusion coefficient',['figuresSptVsWindows/diffCoefModes_band' num2str(iBand) '_good.fig']);
%     
%     %particle density
%     plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.spDensity,...
%         'Particle density',['figuresSptVsWindows/density_band' num2str(iBand) '_good.fig']);
%     
%     %frame-to-frame displacement
%     plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.f2fDispMag2D,...
%     'Frame-to-frame displacement',['figuresSptVsWindows/f2fDisp_band' num2str(iBand) '_good.fig']);
%     
% end
% % close all
% 

