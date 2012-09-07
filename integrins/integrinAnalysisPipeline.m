
%% General

% 1.1 unzip single molecule tif series

% 1.2 convert edge stack to a tif series

% 1.3 remove cell fill frames from single molecule series

%% Single molecule Part 1

% 2.1 crop background area out of single molecule series
%     use: fsmCenter

% 2.2 detect single molecules
%     use: scriptDetectGeneral

% 2.3 track single molecules
%     use: scriptTrackGeneral

% 2.4 run diffusion analysis
%     also put tracks and diffusion analysis in big structures
%     use: scriptCombineTracksDiffAnalysis

%% Cell edge

% 3.1 enhance cell fill using single molecule data
%     use: enhanceEdgeWithSingleMolSignal(?,?,?,400,20)

% 3.2 segment and divide cell mask into windows
%     use: newWindowingCommandsFromHunter

% 3.3 get edge motion characteristics
%     use: classifyEdgeMotion(protSamples,-1,?,?,1)

% 3.4 add movie to list of movies in analysis/120907_MovieInfoStructures

%% Single molecule Part 2

% 4.1 retain tracks >= 5 frames and in cell mask
%     use: scriptRetainTracksCellMask

% 4.2 run diffusion mode decomposition
%     use: scriptDiffModeAnalysis34
%     for this add movie to analysisStruct's in analysis/120702_diffusionModeAnalysis
%     then use: estimDiffModeDividers
%     store results in analysis/120702_diffusionModeAnalysis
%     also generate diffusion mode dividers for classification, as stored
%     in analysis/120702_diffusionModeAnalysis/figuresResults120830/
%     if new movies are of a molecule imaged and analyzed previously, then
%     run diffusion mode decomposition to verify that it has the same modes
%     as the other movies of that molecule, then use the existing diffusion
%     mode dividers.

% 4.3 run diffusion mode classification
%     use: scriptDiffModeClassification

%% Single molecule Part 3

% 5.1 assign tracks to windows
%     use: scriptAssignTracksWindows which calls assignTracks2Windows

% 5.2 calculate "direct track parameters"
%     use: scriptDirectTrackParam which calls trackMotionCharProtrusion

% 5.3 collect particle numbers in small windows
%     use: scriptAssignNumbersWindows which calls assignNumbers2Windows

% 5.4 collect single molecule behavior relative to cell edge activity
%     use: scriptAnalysisSptVsWindows which calls sptRelToActivityOnsetAdaptiveWindows

% 5.5 plot single molecule behavior relative to cell edge activity
%     use: scriptMakePlots which calls plotSptRelToActivityOnsetAdaptiveWindows



