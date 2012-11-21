
%% 1 General

% 1.1 unzip single molecule tif series

% 1.2 convert edge stack to a tif series

% 1.3 remove cell fill frames from single molecule series
%     use: removeSubsetOfFiles

%% 2 Single molecule Part A

% 2.1 crop background area out of single molecule series
%     use: fsmCenter

% 2.2 detect single molecules
%     use: scriptDetectGeneral

% 2.3 track single molecules
%     use: scriptTrackGeneral

% 2.4 run diffusion analysis
%     also put tracks and diffusion analysis in big structures
%     use: scriptCombineTracksDiffAnalysis

%% 3 Cell edge

% 3.1 enhance cell fill using single molecule data
%     use: enhanceEdgeWithSingleMolSignal(?,?,?,400,20)

% 3.2 segment and divide cell mask into windows
%     use: newWindowingCommandsFromHunter

% 3.3 get edge motion characteristics
%     use: classifyEdgeMotion(protSamples,-1,?,?,1)

% 3.4 add movie to list of movies in analysis/120907_MovieInfoStructures

%% 4 Single molecule Part B

% 4.1 retain tracks >= 5 frames and in cell mask
%     use: scriptRetainTracksCellMask

% 4.2 run diffusion mode decomposition
%     if condition is analyzed for 1st time, assemble its movieStruct as
%     stored in analysis/120907_MovieInfoStructures.
%     if there is an existing movieStruct for condition, add movie(s) to
%     it.
%     then use: scriptDiffModeAnalysis34 which stores results in individual
%     directories and collectively in
%     analysis/120702_diffusionModeAnalysis.
%     then use scriptGetDiffModeDividers
%     also generate diffusion mode dividers for classification
%     if new movies are of a molecule imaged and analyzed previously, then
%     run diffusion mode decomposition to verify that it has the same modes
%     as the other movies of that molecule, then use the existing diffusion
%     mode dividers.

% 4.3 run diffusion mode classification
%     use: scriptDiffModeClassification which calls trackDiffModeAnalysis

%% 5 Single molecule Part C

% 5.1 assign tracks to windows
%     use: scriptAssignTracksWindows which calls assignTracks2Windows

% 5.2 calculate "direct track parameters"
%     use: scriptDirectTrackParam which calls trackMotionCharProtrusion

% 5.3 collect particle numbers in small windows
%     use: scriptAssignNumbersWindows which calls assignNumbers2Windows

% 5.4 measure single molecule behavior relative to cell edge activity
%     use: scriptAnalysisSptVsWindows which calls sptRelToActivityOnsetAdaptiveWindows

% 5.5 plot single molecule behavior relative to cell edge activity
%     use: scriptMakePlots which calls plotSptRelToActivityOnsetAdaptiveWindows

%% 6 Single molecule Part D

% 6.1 combine single cell measurements for all cells of same condition
%     use: scriptSptVsWindowsCombineCells & sptRelToActivityMultipleCells

% 6.2 plot collective cell behavior
%     use: scriptMakeCombPlots which calls plotSptRelToActivityOnsetAdaptiveWindows

