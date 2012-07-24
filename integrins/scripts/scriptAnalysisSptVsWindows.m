clear all
close all

% delete edgeMotionChar/*
% delete figuresSptVsWindows/*
% delete particleBehaviorRelActivity.mat

% windowRange = (39:65)';
windowRange = [];
frameRange = [];

% load ../../analysisCellEdgeMod/protrusion_samples/protrusion_samples.mat
% [windowMotionType,windowMotionChar] = classifyEdgeMotion(protSamples,-1,windowRange,[],1,'edgeMotionChar');
% save('edgeMotionChar/edgeMotionChar','windowMotionType','windowMotionChar')
% % close all

load tracksDiffusionLength5InMask.mat
load directTrackChar.mat
load windowsActivityTracks.mat
load diffusionModeClassification.mat
load edgeMotionChar/edgeMotionChar.mat

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[sptPropInWindow,~,~,analysisParam] = ...
    particleBehaviorRelToActivityOnset(tracksFinal,diffAnalysisRes,...
    diffModeAnalysisRes,trackChar,windowsAll,1:400:maxFrame+1,protSamples,...
    5,[1 1; 2 2; 3 3],windowRange,frameRange,windowTrackAssign);

save('particleBehaviorRelActivityAll','sptPropInWindow','analysisParam');

for iBand = 1 : 2
    
    %fraction in each diffusion mode
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.fracMode,...
        'Fraction in diffusion mode',['figuresSptVsWindows/fracModes_band' num2str(iBand) '_all.fig']);
    
    %overall diffusion coefficient
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.diffCoefModeAll,...
        'Overall diffusion coefficient',['figuresSptVsWindows/diffCoefModes_band' num2str(iBand) '_all.fig']);
    
    %particle density
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.spDensity,...
        'Particle density',['figuresSptVsWindows/density_band' num2str(iBand) '_all.fig']);
    
    %frame-to-frame displacement
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.f2fDispMag2D,...
        'Frame-to-frame displacement',['figuresSptVsWindows/f2fDisp_band' num2str(iBand) '_all.fig']);
    
end
% close all

windowRange = windowMotionChar.indxGoodWindow;
frameRange = [];

[sptPropInWindow,~,~,analysisParam] = ...
    particleBehaviorRelToActivityOnset(tracksFinal,diffAnalysisRes,...
    diffModeAnalysisRes,trackChar,windowsAll,1:400:maxFrame+1,protSamples,...
    5,[1 1; 2 2; 3 3],windowRange,frameRange,windowTrackAssign);

save('particleBehaviorRelActivityGood','sptPropInWindow','analysisParam');

for iBand = 1 : 2
    
    %fraction in each diffusion mode
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.fracMode,...
        'Fraction in diffusion mode',['figuresSptVsWindows/fracModes_band' num2str(iBand) '_good.fig']);
    
    %overall diffusion coefficient
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).diffModeAnalysis.diffCoefModeAll,...
        'Overall diffusion coefficient',['figuresSptVsWindows/diffCoefModes_band' num2str(iBand) '_good.fig']);
    
    %particle density
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.spDensity,...
        'Particle density',['figuresSptVsWindows/density_band' num2str(iBand) '_good.fig']);
    
    %frame-to-frame displacement
    plotParticleBehaviorRelToActivityOnset(sptPropInWindow(iBand).directAll.f2fDispMag2D,...
    'Frame-to-frame displacement',['figuresSptVsWindows/f2fDisp_band' num2str(iBand) '_good.fig']);
    
end
% close all


