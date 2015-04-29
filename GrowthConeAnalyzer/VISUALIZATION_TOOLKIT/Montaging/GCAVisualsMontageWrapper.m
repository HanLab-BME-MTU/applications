function [ output_args ] = GCAVisualsMontageWrapper( projList )
% Wrapper function that takes in a cell array of movieData (project paths) 
% Checks to see if the raw image is run - if not runs it
% Checks to see if the overlay is run - if not runs it 
% Checks for the montage directory - if not present runs montaging and
% makes a movie. 
% 
% 
%
%
% Test for which projects have already been run
idxNeedToRunRaw = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'Raw'])==0,projList);
idxNeedToRunOverlay = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'Overlaysmono'])==0,projList);

% find what needs to be run
list{1} = projList(idxNeedToRunRaw);
list{2} = projList(idxNeedToRunOverlay);

visFuncH{1} = 'GCAVisualsInvertedRawMovie';
visFuncH{2} = 'GCAVisualsMakeOverlaysFilopodiaOldInputMovie';

% for all the movies where the overlay directories that you want have not
% be run : run them now!
success = cellfun(@(x,y) GCAVisualsWrapper(x,y),list,visFuncH);




% check which projects do not have project directories
idxNeedToRun = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'montage'])==0,projList);
cellfun(@(x) display(['Montage Directory found for ' x ' Skipping']),projList(~idxNeedToRun));
projListRunMontage = projList(idxNeedToRun);


for iProj = 1:numel(projListRunMontage)
    load([projList{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']);
    
    % make the inputDirectories for the montages.
    inputDirs{1} = [MD.outputDirectory_ filesep ...
         'Visualization_Overlays' filesep 'Raw'];
    inputDirs{2} = [MD.outputDirectory_ filesep ...
         'Visualization_Overlays' filesep 'Overlaysmono'];
    
    % make the montage directory
    outputDir = [MD.outputDirectory_ filesep ...
         'Visualization_Overlays' filesep 'montage'];
     if ~isdir(outputDir) 
         mkdir(outputDir) 
     end 
    
    if size(projList,2)>1;
        % get the name
        nameMovie = projList{iProj,2};
    else
        nameMovie =[];
    end
    GCAVisualsMontagingMovie(MD,inputDirs,outputDir,nameMovie);
end
end



