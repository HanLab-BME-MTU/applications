function [ collectedData ] = plusTipScreeningFunction(varargin)

% a plusTipTracker AD-HOC function for testing for parameter hit
% consistency among multiple experiments. (used for the VanVactors screening paper)
% note that in the future
% dataset array generation and manipulation is likely a more efficient
% set-up for data manipulation of this mag- though hind-site is 20/20. :(
%
% -MB

% Input:
% createDayList : 1 if you would like to create the dayList. 0 it will ask
% you to load. Default = 0
%
% ID: nx1 cell array of IDs of conditions for which you would like to find
% hits (sorry there is no easy way to make this however, once you do you
% can save it and load the data) Default  = [] it will ask you to load a
% file
%
% analysisFolder: string that tells which analysis folder under each day
% you would like to collect if you have run the analysis through the
% . ex 'analysis_01_29_PooledReclass' (default = 'Analysis')
%
% stringency = the p-value cut-off for which you would like to collect
%
% percentDays = the percentage of days that need to be significant in order
% to be recorded. (default = 0.5)
%
% makeTextFiles = set to 1 to make text (default = 1)
%
% makeCollectedPlots = set to 1 to make these plots. (default = 0)
%
%
%% Check Input
ip = inputParser;
ip.addParamValue('createDayList',0,@isscalar);
ip.addParamValue('ID',@(x)iscell(x) ||isempty(x));
ip.addParamValue('analysisFolderName','Analysis',@ischar);
ip.addParamValue('stringency',0.05,@isscalar);
ip.addParamValue('percentDays',0.5,@isscalar);
ip.addParamValue('makeTextFiles',1,@isscalar);
ip.addParamValue('makeCollectedPlots',0,@isscalar);

ip.parse(varargin{:});
createDayList = ip.Results.createDayList;
ID = ip.Results.ID;
analysisFolderName = ip.Results.analysisFolderName;
stringency = ip.Results.stringency;
percentDays = ip.Results.percentDays;
makeTextFiles = ip.Results.makeTextFiles;
makeCollectedPlots = ip.Results.makeCollectedPlots;


if isempty(ID)
    IDFilename  =  uigetfile('*.mat','Please Select a File Containing the ID Cell');
    s = load(IDFilename);
    ID = s.ID;
end
% Create OR LOAD DAY LIST

if createDayList == 1
    % make a List of Day Paths
    flag = 0;
    dayFileCount = 1;
    while flag == 0
        dayPath =uigetdir(pwd,'Select Directory Under Which Screen Stat File is Located');
        
        % save to day path list
        dayPathList{dayFileCount,1} = dayPath;
        
        choice = questdlg('Select Another Day?','Select Another Day','Yes','No','No');
        
        switch choice
            case 'Yes'
                flag = 0 ;
                dayFileCount = dayFileCount + 1;
            case 'No'
                flag = 1;
        end
    end
    
    [saveDir ]  = uigetdir(pwd,'Where would you like to store the day list');
    temp=inputdlg({'Enter file name:'},'',1);
    save([saveDir filesep temp{1}],'dayPathList');
    
else
    [file2Load, path] = uigetfile(pwd,'Select Day List to Load');
    temp = load([path filesep file2Load]);
    dayPathList = temp.dayPathList;
    
    
end % end if Create Day List

[saveDir] = uigetdir(pwd,'Where would you like to store the output?');
% stat file hard input
nameOfStatFile1 = [ analysisFolderName filesep 'perCell' ];
nameOfStatFile2 = analysisFolderName;

%%
collectedData = struct();
[days ~ ] = size(dayPathList);

for iDay = 1:days
    % load hit structure %%% NOTE THIS IS WHERE YOU ARE CHOOSING WHICH MAT
    % TO LOAD Eventually have as an input
    if isdir([dayPathList{iDay,1} filesep nameOfStatFile1])
        file2load =[dayPathList{iDay,1} filesep nameOfStatFile1 filesep 'hitsTest2.mat']; % here the hitsTest2 should correspond to the perm t-test
    else
        file2load = [dayPathList{iDay,1} filesep nameOfStatFile2 filesep 'hitsTest2.mat'];
    end
    
    % New work-up (for Mijung)
    %file2load = [dayPathList{iDay,1} filesep 'perCell' filesep 'hitsTest2.mat'];
    
    
    temp = load(file2load);
    
    
    
    
    %% Extract Data for all Parameters and cluster
    % grab the title to carry over.
    
    %if iDay == 1
    %   title = data(1,:);
    %end
    % get the param names
    if isfield(temp,'hits2') ==1
        statNames = fieldnames(temp.hits2);
    else
        statNames = fieldnames(temp.hits);
    end
    % not the most elegant way to do (maybe modify in the future)
    if iDay == 1
        statsFinal = statNames;
    else
        statsFinal = [statsFinal;setdiff(statsFinal,statNames)];
    end
    totStats = size(statNames,1);
    totIDs = length(ID);
    % extract data for each ID
    
    for iStat= 1:totStats
        
        
        for iID = 1:totIDs
            % find the info that you want
            if isfield(temp ,'hits');
                data = temp.hits.(char(statNames(iStat)));
                
            else
                if isfield(temp.hits2,statNames(iStat))==1;
                    data = temp.hits2.(char((statNames(iStat))));
                end
            end
            
            
            
            
            % take out grp
            dataNames = cellfun(@(x) strrep(x,'grp_',''),data(:,1),'uniformoutput',0);
            % idx of all
            
            %using sting compare here fairly sensitive maybe better way to do this
            idx = cell2mat(cellfun(@(x) strncmpi(x, ID{iID,1},length(ID{iID,1})) ,dataNames(:,1),'uniformoutput',0));
            
            % just get the info you want save for each (overwrite for each)
            
            data2Save = data(idx,:);
            if ~isempty(data2Save)
                
                data2Save = horzcat(data2Save,dayPathList{iDay}); % save the directory so can use this later to plot hits
                
                % compile each day (not a good way to do this change in future)
                
                if ( ~isfield(collectedData,(char(statNames(iStat))))== 1 || ~isfield(collectedData.(char(statNames(iStat))),(char(ID{iID,1})) )== 1);
                    
                    collectedData.(char(statNames(iStat))).(char(ID{iID})) = data2Save;
                else
                    
                    collectedData.(char(statNames(iStat))).(char(ID{iID}))= [collectedData.(char(statNames(iStat))).(char(ID{iID,1})) ; data2Save];
                    
                end % end if isfield
            end
            %
            
            
            
            
        end  % stats
        
    end % end iID
end % for iDay
collectedDataFiltered = struct();
collectedDataFilterFinal = struct();
collectedDataFilterAll = struct();

% filter collectedData
for iStat = 1:length(statsFinal)
    for iID = 1:size(ID,1);
        
        data  = collectedData.(char(statsFinal(iStat))).(char(ID{iID,1}));
        [total ~ ] = size(data);
        hits = data((cell2mat(data(:,4)) < stringency),:);
        
        
        sig = size(hits,1) ;
        if (sig/total >= percentDays)
            collectedDataFilterAll.(char(statsFinal(iStat))).(char(ID{iID,1})) = data;
            nonhits = data((cell2mat(data(:,4))>stringency),:);
            titleHits = cell(1,5);
            titleHits{1,1} = 'Hit Days';
            %                 titleHits{2,1} = 'Name';
            %                 titleHits{2,2} = 'Value';
            %                 titleHits{2,3} = 'Normalized by Control';
            %                 titleHits{2,4} = 'P-value';
            %                 titleHits{2,5} = 'Path';
            
            
            values = cell2mat(hits(:,3));
            percentChange = (values-1)*100;
            
            testForSign = sign(percentChange);
            if abs(sum(testForSign))<length(percentChange) ;
                % check to see which is the problem
                badIdx = find(testForSign~=mode(testForSign)) ;
                bad = hits(badIdx,:)  ;
                nonhits = [nonhits;cell(1,5);bad] ;
                % bad way to do this but to save time
                for i = 1:length(badIdx)
                    hits{badIdx(i),1} = 'Warning Significant Hit In Opposite Direction';
                    hits{badIdx(i),2} = [];
                    hits{badIdx(i),3} = [];
                    hits{badIdx(i) ,4} = [];
                end
                
            end
            diff = size(hits,1) - size(nonhits,1) ;
            
            add = cell(diff+7,5);
            
            titleNonHits = titleHits;
            titleNonHits{1,1} = 'Non Hit Days';
            
            
            
            nextLine = size(hits,1)+1;
            hits{nextLine,1} = 'Average' ;
            hits{nextLine,2} = mean(cell2mat(hits(:,2)));
            hits{nextLine,3} = mean(cell2mat(hits(:,3)));
            hits{nextLine,4} = mean(cell2mat(hits(:,4)));
            hits{nextLine+1,1} = 'Average Percent Change Excluding Hits';
            values = cell2mat(hits(:,3));
            percentChangeHitValues= (values-1)*100;
            
            hits{nextLine+1,3} = mean(percentChangeHitValues);
            hits{nextLine+2,1} = 'SEM Percent Change Excluding NonHit Days';
            hits{nextLine+2,3} = std(percentChangeHitValues)/sqrt(length(percentChangeHitValues));
            
            nonHitValues = cell2mat(nonhits(:,3));
            percentChangeNonHitValues = (nonHitValues-1)*100;
            percentChangeAllValues = [percentChangeNonHitValues;percentChangeHitValues];
            hits{nextLine+4,1} = 'Average Percent Change Including Non Hit Days';
            hits{nextLine+4,3} = mean(percentChangeAllValues);
            hits{nextLine+5,1} = 'SEM Percent Change Including NonHit Days';
            hits{nextLine+5,3} = std(percentChangeAllValues)/sqrt(length(percentChangeAllValues));
            hits{nextLine+6,1} = []; % make extra line
            
            hits = [titleHits;hits];
            nonhits = [titleNonHits; nonhits;add];
            hits  = [hits nonhits];
            collectedDataFilterFinal.(char(statsFinal(iStat))).(char(ID{iID,1})) = hits;
            
            
        else
            % don't record data
        end
        
    end
end




save([saveDir filesep 'collectedData'],'collectedData');
save([saveDir filesep 'collectedDataFilterAll'],'collectedDataFilterAll');
save([saveDir filesep 'collectedDataFilterFinal'],'collectedDataFilterFinal');


% params = fieldnames(collectedDataFiltered) ;
% % this is probably still set up poorly but it works
% for iParam = 1:numel(params)
%
%     hitsForParam = fieldnames(collectedDataFiltered.char(params(iParam)));
%     for iHit = 1:length(hitsForParam)
%       collectHits = collectedData.(char(params(iParam)).char(hitsForParam(iHit)).hitsForParam(iHit);
%
%     end
%
% end



%% save text file and plots of data

% make independent folder with hit information given that parameter
if makeTextFiles== 1;
    
    if ~isdir([saveDir filesep 'Hit_Information'])
        
        mkdir([saveDir filesep 'Hit_Information'])
    end
    
    %     if ~ isdir([saveDir filesep 'CollectedTextFiles']);
    %     mkdir([saveDir filesep 'CollectedTextFiles']);
    %     end
    %
    %for each stat that was a hit
    
    for iStat = 1:size(fieldnames(collectedDataFilterAll),1)
        %make txt file of hits for that param
        
        statsFiltered = fieldnames(collectedDataFilterAll);
        statName = char(statsFiltered(iStat));
        statsFile = [saveDir filesep 'Hit_Information' filesep statName];
        
        if ~isdir(statsFile)
            mkdir(statsFile);
        end
        
        
        output = struct2cell(collectedDataFilterFinal.(char(statsFiltered(iStat))));
        output  = vertcat(output{:});
        
        
        % statsFile = ['hits_' char(statsFiltered(iStat)) ] ;
        %title{1,1} = 'Hit Name';
        %title{1,2} = 'Mean Value of Cellular Population';
        %title{1,3} = 'Mean Value of Cellular Population Normalized to Control';
        %title{1,4} = 'P-value';
        
        fid=fopen([statsFile filesep statName 'txt'],'w+');
        
        %fprintf(fid,'\t%s',char(title));
        nameOfParam = char(statsFiltered(iStat));
        
        % condNames = output(:,1);
        % values = output(:,2:4);
        sepCol = cellfun(@(x)find(strcmp(x,output(1,:))),{'Hit Days','Non Hit Days'});
        sepCol(end+1)=size(output,2)+1;
        
        fprintf(fid,'\n%s\n',nameOfParam);
        for i=1:size(output,1)
            for j =1:numel(sepCol)-1
                condNames = output(i,sepCol(j));
                values = output(i,sepCol(j)+1:sepCol(j+1)-2);
                paths = output(i,sepCol(j+1)-1);
                
                
                fprintf(fid,'%s\t',condNames{:});
                fprintf(fid,'%g\t',values{:});
                fprintf(fid,'%s\t',paths{:});
                fprintf(fid,'\t');
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        %%
        if makeCollectedPlots == 1
            hitnames = fieldnames(collectedDataFilterAll.(char(statsFiltered(iStat))));
            % make collected plots
            
            % load all the days that were hits
            
            
            dataStruct = struct();
            
            for iCond = 1:length(hitnames)
                charName = strrep(char(hitnames{iCond}),'_','');
                
                saveDirStatCond = [statsFile filesep charName];
                if ~isdir(saveDirStatCond)
                    mkdir(saveDirStatCond)
                end
                
                
                data =  collectedDataFilterAll.(char(statsFiltered(iStat))).(char(hitnames{iCond}));
                count = 0;
                for iDay = 1:length(data(:,1))
                    
                    
                    groupList = load([char(data(iDay,5)) filesep 'groupList.mat']);
                    groupList = groupList.groupList;
                    
                    % assumes control is named such should change to be first group?
                    groupListControl = groupList(cell2mat(cellfun(@(x) strncmpi(x, 'Cont',3) ,groupList(:,1),'uniformoutput',0)),:);
                    % condNames =  cellfun(@(x) strrep(x,'grp_',''),char(hitnames{iCond}),'uniformoutput',0);
                    groupListHit = groupList(cell2mat(cellfun(@(x) strncmpi(x,char(hitnames{iCond}),3),groupList(:,1),'uniformoutput',0)),:);
                    
                    
                    
                    groupList = [groupListControl;groupListHit];
                    
                    
                    % in the new workflow will have saved this already: however here just
                    % do again
                    groupData = plusTipExtractGroupData(groupList,1);
                    
                    
                    
                    
                    
                    cellValuesCont = ...
                        cellfun(@(x) x.(char(statsFiltered{iStat})),groupData.stats{1});
                    cellValuesCont = abs(cellValuesCont);
                    group  = (2*count+1)*ones(1,length(cellValuesCont));
                    scatter(group,cellValuesCont,100,'b','Filled');
                    hold all
                    plot(2*count+1,mean(cellValuesCont),'+','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10) ;
                    dataStruct(iDay).cellValuesCont = cellValuesCont';
                    
                    cellValuesHit = ...
                        cellfun(@(x) x.(char(statsFiltered{iStat})),groupData.stats{2});
                    cellValuesHit = abs(cellValuesHit);
                    
                    group = (2*count+2)*ones(1,length(cellValuesHit));
                    scatter(group,cellValuesHit,100,'r','Filled');
                    plot(2*count+2,mean(cellValuesHit),'+','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10),
                    dataStruct(iDay).cellValuesHit = cellValuesHit';
                    count = count+3;
                    
                    % make movieList
                    % find(cellValuesCont == max(cellValuesCont))
                    
                end
                
                nameStat = strrep(char(statsFiltered{iStat}),'_',' ');
                ylabel(nameStat,'fontsize',20);
                
                maxX = (2*count+3);
                minVal =  min([vertcat(dataStruct(:).cellValuesHit);vertcat(dataStruct(:).cellValuesCont)]) ;
                maxVal =  max([vertcat(dataStruct(:).cellValuesHit);vertcat(dataStruct(:).cellValuesCont)]);
                diff = maxVal - minVal;
                
                minY = minVal-0.2*(diff);
                maxY = maxVal+0.2*(diff);
                axis([ 0 maxX minY maxY ]);
                
                
                
                legend('Control','Control Mean', charName,[ charName ,'Mean']);
                saveas(gcf,[saveDirStatCond filesep 'cellDistCollected_' char(hitnames{iCond}) '.fig']);
                saveas(gcf,[saveDirStatCond filesep 'cellDistCollected_' char(hitnames{iCond}) '.eps'],'psc2');
                hold off
                
                close(gcf)
                
            end % iCond
            
        end
        
    end  % iStat
    
end

