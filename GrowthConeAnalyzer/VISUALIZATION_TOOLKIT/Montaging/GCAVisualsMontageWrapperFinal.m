function [ output_args ] = GCAVisualsMontageWrapperFinal( projList, visType ,collect)
% Wrapper function that takes in a cell array of movieData (project paths)
% Checks to see if the raw image is run - if not runs it
% Checks to see if the overlay is run - if not runs it
% Checks for the montage directory - if not present runs montaging and
% makes a movie.
% visType: Cell array of types of visualizations to run (Required)
% collect.on = true  ( param: optional)
% collect.path = name of path to save (if no path ask the user where to
% save)
% collect.type = default(montage)
% montage.on:  make a montaging directory
% montage.visualsToMontage:
%

% default visualization type and functions
if nargin<2 || isempty(visType)
    visType{1} =  'Raw';
    visType{2} =  'Overlaysmono';
    visFuncH{1} = 'GCAVisualsInvertedRawMovie';
    visFuncH{2} = 'GCAVisualsMakeOverlaysFilopodiaOldInputMovie';
end

montage.on = true;
% for all the visuals
for iVisual = 1:numel(visType)
    % Test for which projects have already been run
    % note will eventually replace by checking the movieData after you write
    % processes for all of these but for now need to just check for the
    % existence of file...
    
    
    % test all the project for given visualization
    idxNeedToRunC = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep visType{iVisual}])==0,projList);
    % get the list of projects that still need to be run
    listC = projList(idxNeedToRunC);
    % initial success flag for entire project list based on the idxNeedToRun
    success = ~idxNeedToRunC; % flag where ok
    indexNeedToRunC=  find(idxNeedToRunC); % i know this is stupid but find the original index
    
    
    % if the list isn't empty run them  now and see if they are successful
    if ~isempty(listC)
        afterFix = GCAVisualsWrapper(listC,visFuncH{iVisual}); % will be lenght of projs
        for i = 1:length(afterFix)
            success(indexNeedToRunC(i)) = afterFix(i);  % update the success flag based on the success of the run
        end
    end
    successFinal{iVisual} = success;% % record the success final for each visual should be the length of the original projList
    
end % iRun (type of visual

% get the final success rate



% % for all the movies where the overlay directories that you want have not
% % be run : run them now! success will be a cell array of the lists for each
% % visualization that needed to be run
% success = cellfun(@(x,y) GCAVisualsWrapper(x,y),list,visFuncH,'uniformoutput',0); % will wrap through all the project list
% find all the projects that were not a success and remove them.


% associated with all the visualization handles.
%% filter project for montaging by success
if montage.on == true
    %readyToGo{iRun};
    %y = list{iRun} ;
    %% Filter by if the files exist for the montage
    
    successAllVisuals = horzcat(successFinal{:});
    % test to see which projects you can actually make a full montage given
    % the ability to make the overlay requested.- ie did not get a success
    % flag for all visualization types.
    idxProjectNoMont =  sum(successAllVisuals,2)./length(successAllVisuals(:,2)) ~= 1;
    projListNo = projList(idxProjectNoMont,:);
    projListForMontaging = projList(~idxProjectNoMont,:);
    cellfun(@(x) display(['Cannot make montage of ' x]),projListNo);
    % if collectDir ==1
    %% Now that you have all the files you need start montaging (eventually make input)
    
    
    % check which projects do not have project directories
    
    
    idxNeedToRun = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'montage'])==0,projListForMontaging(:,1));
    cellfun(@(x) display(['Montage Directory found for ' x ' Skipping']),projListForMontaging(~idxNeedToRun,1));
    projListRunMontage = projListForMontaging(idxNeedToRun,:);
    
    
    for iProj = 1:numel(projListRunMontage)
        load([projListRunMontage{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']);
        
        
        inputDirs =cellfun(@(x) [ MD.outputDirectory_ filesep 'Visualization_Overlays' filesep x],visType,'uniformoutput',0);
        
        % make the inputDirectories for the montages.
        %     inputDirs{1} = [MD.outputDirectory_ filesep ...
        %          'Visualization_Overlays' filesep 'Raw'];
        %     inputDirs{2} = [MD.outputDirectory_ filesep ...
        %          'Visualization_Overlays' filesep 'Overlaysmono'];
        
        % make the montage directory
        outputDir = [MD.outputDirectory_ filesep ...
            'Visualization_Overlays' filesep 'montage'];
        if ~isdir(outputDir)
            mkdir(outputDir)
        end
        
        if size(projListRunMontage,2)>1;
            % get the name
            nameMovie = projListRunMontage{iProj,2};
        else
            nameMovie =[];
        end
        GCAVisualsMontagingMovie(MD,inputDirs,outputDir,nameMovie);
        
        
    end % iProj
    %% Collect the files if the user desires.
    for iProj = 1:numel(projListForMontaging(:,1))
        
        if collect.on == true
            
            % find the mp4 file and move it to the collect dir
            source = searchFiles('.mp4',[],[projListForMontaging{iProj} filesep 'ANALYSIS' filesep ...
                'Visualization_Overlays' filesep 'montage'],0);
            
            
            dest = cellfun(@(x)[ collect.path filesep x],source(:,1),'uniformoutput',0);
            
            cellfun(@(x) display(['Collected Montage for Project ' x]),projListForMontaging);
            arrayfun(@(x) movefile([source{x,2} filesep source{x,1}],dest{x}),1:numel(source(:,1)));
            %  GCAVisualse
        end % if collect.on
    end % iProj
    %% If can't make montage move over just the raw image.
    for iProj = 1:numel(projListNo)
        cDir = [projListNo{iProj,1} filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'Raw'];
        % search for .mp4s
        source = searchFiles('.mp4',[],cDir,0);
        load([projList{iProj,1} filesep 'movieData.mat']); 
        if size(projListNo,2) > 1 % you have movie Names
            movieName = projListNo(iProj,2);
        else
            movieName = 'movie';
        end
        if isempty(source) % didn't find a movie make the movie but need to crop first 
            cropImageToMakeMovieOnWindows(MD); 
            
            
            movefile([cropDir filesep movieName '.mp4'],[collect.path filesep movieName '.mp4']);
        else
            movefile([cDir filesep source{1,:} '.mp4'],[collect.path filesep movieName '.mp4']);
        end % isempty(source)
        % search Again
        
        
        
        
    end % iProj
    
    
    
end % montage.on == true




