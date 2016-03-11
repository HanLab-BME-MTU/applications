function [ output_args ] = GCAAnalysisToolsMakeToPlotFileFromGroups(searchDirectory,screenByNotes,findExcludeFolders)
%% 
%  searchDirectory: the directory with the larger knockdowns to choose from
%  
% 
% screenByNotes: if true will look in the notes directory and see if you 
% selected it for analysis 
% findExcludeFolder: currently set-up filters those folders which are not 
% just numbers in case there is extra crap in there- here I include those. 
% 

excludeReport = []; 
% get everything with KD
% fitlerByNotes: will check the notes flag ...
s = dir(searchDirectory);
idxKeep =  arrayfun(@(x) ~isempty(regexp(x.name,'KD','ONCE')),s);
s = s(idxKeep);
groupName = arrayfun(@(x) [x.name],s,'uniformoutput',0);

idxInclude  = listSelectGUI(groupName,[],'move');

groupName = groupName(idxInclude) ;

nGroups = numel(groupName);
% for each group
for iGroup = 1:nGroups
    toPlot.info.names{iGroup} = groupName{iGroup};
    % for each group look for experimental day folders
    directoriesExpDay = dir([searchDirectory filesep groupName{iGroup}]);
    idxDir = vertcat(directoriesExpDay(:).isdir);
    directoriesExpDay = directoriesExpDay(idxDir);
    % keep only those experimental day folders with numbers
   
    idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesExpDay);
    directoriesExpDay = directoriesExpDay(idxKeep);
    
    %directoriesExpDay = directoriesExpDay(idxKeep);
    directoriesExpDay =   arrayfun(@(x) [searchDirectory filesep groupName{iGroup} filesep x.name],directoriesExpDay,'uniformoutput',0);
    % for each day folder look for movie folders (a flag will be in the
    % notes on whether not to currently include or not)
    for iDay = 1:length(directoriesExpDay);
        directoriesMovie = dir(directoriesExpDay{iDay});
        idxDir = vertcat(directoriesMovie(:).isdir);
        directoriesMovie = directoriesMovie(idxDir);
        % keep only the number folders in case I put crap in there.
        idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesMovie);
     
        
        if findExcludeFolders == true
        % Find EXCLUDE 
        idxKeep2 = arrayfun(@(x) ~isempty(regexpi(x.name,'exclude')),directoriesMovie); 
        idxKeep = (idxKeep|idxKeep2); 
        end 
        
        directoriesMovie = directoriesMovie(idxKeep);
        
        forProjList  = arrayfun(@(x) [directoriesExpDay{iDay} filesep  x.name],directoriesMovie,'uniformoutput',0);
        if screenByNotes == true
        % test the notes to see if should include
        notesAll = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'notes.mat']), forProjList);
        idxInclude = arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,'yes')),1:length(notesAll));
        
        excludedFiles = forProjList(~idxInclude); 
        if ~isempty(excludedFiles); 
         reasons  = arrayfun(@(x) notesAll(x).notes.ExcludeAtStep,1:length(notesAll),'uniformoutput',0);
         reasons =  reasons(~idxInclude); 
         excludeReport{iDay} = [excludedFiles reasons']; 
        end 
        
        else 
            idxInclude = true(size(forProjList,1),1); % include them all 
        end 
        forProjList = forProjList(idxInclude);
        namesProj = cell(numel(forProjList),1);
        
        for iProj = 1:numel(forProjList)
            
            [~,date] = upDirectory(forProjList{iProj},2,1);
            [~,num] = upDirectory(forProjList{iProj},1,1);
            [~,name] = upDirectory(forProjList{iProj},3,1);
            namesProj{iProj,1}= [name '_' date '_' num];
        end
        projListC = [forProjList namesProj];
        projListGroup{iDay} = projListC;
    end
    
    
    toPlot.info.projList{iGroup} = vertcat(projListGroup{:});
    nProjs = size(toPlot.info.projList{iGroup},1);
    grouping{iGroup} = iGroup.*ones(nProjs,1);
  clear projListGroup
end
toPlot.info.grouping = vertcat(grouping{:});


save('toPlotGroup.mat','toPlot');
save('ExcludeReport.mat','excludeReport'); 

end

