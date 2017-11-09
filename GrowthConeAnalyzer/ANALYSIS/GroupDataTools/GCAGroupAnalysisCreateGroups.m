function [toPlot] = GCAGroupAnalysisCreateGroups(searchDirectory,varargin)

%% GCAGroupAnalysisCreateGroups
%
% Note: this file is a bit designed exclusively for my needs... should
% figure out how to revamp before release. 

%  searchDirectory: the directory with the larger knockdowns folders to choose from
%                   I set it up currently so that all the group names
%                   required start with KD as an identifier
%
% screenByNotes: if true will look in the notes directory and see if you
%               selected it for analysis
%
% switched name to GCAGroupAnalysisCreateGroups from GCAAnalysisToolsMakeToPlotFileFromGroupsMod.m on 20151027
%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;
%ip.KeepUnmatched = true;

ip.addRequired('searchDirectory');
ip.addParameter('screenByNotes',true);
ip.addParameter('searchTerm',[]); 
ip.addParameter('OutputDirectory',pwd);
ip.addParameter('mainFig',false); 
ip.addParameter('biosensors',false); 

ip.parse(searchDirectory,varargin{:});

screenByNotes = ip.Results.screenByNotes; 

%% 
excludeReport = [];
% get everything with KD
% fitlerByNotes: will check the notes flag ...
s = dir(searchDirectory);
if ip.Results.biosensors
    idxKeep = arrayfun(@(x) ~isempty(regexp(x.name,'Biosensor','ONCE')),s);
else
    idxKeep =  arrayfun(@(x) ~isempty(regexp(x.name,'KD','ONCE')),s);
end
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
    
    if ip.Results.biosensors
        idxKeep = arrayfun(@(x) ~isempty(regexp(x.name,'KD','ONCE')),directoriesExpDay);
        directoriesExpDay = directoriesExpDay(idxKeep);
        %x =  arrayfun(@(x) x.name,directoriesExpDay,'uniformoutput',0); 
        % for now assume only 1 perturbation condition 
        pertCond = directoriesExpDay.name;  
        clear('directoriesExpDay');  
        directoriesExpDay = dir([searchDirectory filesep groupName{iGroup} filesep pertCond]);   
    end
    
    idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesExpDay);
    directoriesExpDay = directoriesExpDay(idxKeep);
    
    %directoriesExpDay = directoriesExpDay(idxKeep);
    if ip.Results.biosensors
        directoriesExpDay =   arrayfun(@(x) [searchDirectory filesep groupName{iGroup} filesep pertCond filesep x.name],directoriesExpDay,'uniformoutput',0);
    else
        directoriesExpDay =   arrayfun(@(x) [searchDirectory filesep groupName{iGroup} filesep x.name],directoriesExpDay,'uniformoutput',0);
    end
    %
    load([searchDirectory filesep groupName{iGroup} filesep 'groupColor.mat'])
    % for each day folder look for movie folders (a flag will be in the
    % notes on whether not to currently include or not)
    for iDay = 1:length(directoriesExpDay);
        directoriesMovie = dir(directoriesExpDay{iDay});
        idxDir = vertcat(directoriesMovie(:).isdir);
        directoriesMovie = directoriesMovie(idxDir);
        % keep only the number folders in case I put crap in there.
        idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesMovie);
        
        directoriesMovie = directoriesMovie(idxKeep);
        
        forProjList  = arrayfun(@(x) [directoriesExpDay{iDay} filesep  x.name],directoriesMovie,'uniformoutput',0);
        if screenByNotes == true
            % test the notes to see if should include
            notesAll = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'notes.mat']), forProjList);
            %idxInclude = arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,'BranchSensAnal')),1:length(notesAll));
            idxInclude = arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,'yes')),1:length(notesAll));
            
            if ip.Results.mainFig 
                idxInclude = arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,'20151104ReRun')),1:length(notesAll)); 
                %idxInclude = idxInclude & idxInclude2 ; 
            end 
            
            
            if ~isempty(ip.Results.searchTerm)
             idxInclude = arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,ip.Results.searchTerm)),1:length(notesAll)); 
            end 
            excludedFiles = forProjList(~idxInclude);
            if ~isempty(excludedFiles);
                reasons  = arrayfun(@(x) notesAll(x).notes.ExcludeAtStep,1:length(notesAll),'uniformoutput',0);
                reasons =  reasons(~idxInclude);
                excludeReportDay{iDay} = [excludedFiles reasons'];
            else
                excludeReportDay{iDay} = [] ;
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
            if ip.Results.biosensors
                [~,biosensor] = upDirectory(forProjList{iProj},4,1);
                namesProj{iProj,1} = [biosensor '_' name '_' date '_' num];
            else
                namesProj{iProj,1}= [name '_' date '_' num];
            end
        end
        projListC = [forProjList namesProj];
        projListGroup{iDay} = projListC;
    end
    if screenByNotes == true;
        excludeReport{iGroup} = vertcat(excludeReportDay{:});
        clear excludeReportDay
    end
    toPlot.info.projList{iGroup} = vertcat(projListGroup{:});
    nProjs = size(toPlot.info.projList{iGroup},1);
    grouping{iGroup} = iGroup.*ones(nProjs,1);
    toPlot.info.color{iGroup} = groupColor.single;
    toPlot.info.colorShades{iGroup} = groupColor.mult;
    
    clear projListGroup
end
toPlot.info.grouping = vertcat(grouping{:});

if ~isempty(ip.Results.OutputDirectory)
    save([ip.Results.OutputDirectory filesep 'toPlotGroup.mat'],'toPlot');
    save([ip.Results.OutputDirectory filesep 'ExcludeReport.mat'],'excludeReport');
end

