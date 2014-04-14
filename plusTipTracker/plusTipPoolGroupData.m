function plusTipPoolGroupData(groupData,varargin)
% plusTipPoolGroupData pools plus tip data from multiple projects in groups
%
% SYNOPSIS:  [groupData]=plusTipPoolGroupData(groupList,saveDir,doBtw,doWtn,doPlot,remBegEnd)
%
% INPUT:
% groupList : output of plusTipPickGroups, nProj x 2 cell array where the
%             first column contains the group identifier for each group,
%             and the second column contains the project path
% saveDir   : path to output directory
% doBtw     : 1 to pool data from projects in the same group, 0 to skip
% doWtn     : 1 to save text file of distributions from all projects in
%             each group in a "withinGroupComparisons" folder
% doPlot    : 1 to make histograms and boxplots for within and/or between
%             group data
% remBegEnd : 1 to remove tracks existing at the beginning
%             or end of the movie
%
% OUTPUT:
% groupData : structure containing group information and fields for 9
%             distributions: growth speed (gs), fgap speed (fs), bgap speed
%             (bs), growth lifetime (gl), fgap lifetime (fl), bgap lifetime
%             (bl), growth displacement (gd), fgap displacement (fd), and
%             bgap displacement (bd).

% Input check
ip=inputParser;
ip.addRequired('groupData',@(x) isstruct(x) || iscell(x) || isempty(x));
ip.addOptional('saveDir',[],@ischar);
ip.addOptional('doWtn',1,@isscalar);
ip.addOptional('doPlot',1,@isscalar);
ip.parse(groupData,varargin{:});
saveDir=ip.Results.saveDir;
doWtn=ip.Results.doWtn;
doPlot=ip.Results.doPlot;

% Call extract groupData function if empty or string input
if isempty(groupData) || iscell(groupData)
   groupData=plusTipExtractGroupData(groupData); 
end

% Launch interface if empty directory
if isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for pooled group data.');
end

nGroups = numel(groupData.M);
% Within-group comparison
if doWtn
    for iGroup = 1:nGroups
        % Set up output directory
        wtnDir=[saveDir filesep 'withinGroupComparisons' filesep groupData.names{iGroup}];
        if isdir(wtnDir), rmdir(wtnDir,'s'); end
        mkdir(wtnDir);
        
        % write out speed/lifetime/displacement distributions into a text file
        stackedM =  vertcat(groupData.M{iGroup}{:});
        dlmwrite([wtnDir filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd_' ...
            groupData.names{iGroup} '.txt'], stackedM, 'precision', 3,...
            'delimiter', '\t','newline', 'pc');
        
        % Names for each movie in iGroup
        if ~isempty(groupData.movieNames{iGroup})
            wtnGrpNames = groupData.movieNames{iGroup};
        else
            wtnGrpNames = arrayfun(@(x) [groupData.names{iGroup} '_' num2str(x)],...
                1:numel(groupData.dataMat{iGroup}),'Unif',0);
        end

        % Write stats results into a text file
        statsFile = [wtnDir filesep 'Stats.txt'];
        statsNames = fieldnames(groupData.stats{iGroup}{1});
        statsData= cellfun(@struct2cell,groupData.stats{iGroup},'Unif',false);
        statsData =horzcat(statsData{:});
        
        pooledStatsNames = fieldnames(groupData.pooledStats{iGroup});
        pooledStatsData = struct2cell(groupData.pooledStats{iGroup});
        assert(all(ismember(pooledStatsNames, statsNames)));
        
        % Save pooled stats
        fid=fopen(statsFile,'w+');
        fprintf(fid,'\t%s',wtnGrpNames{:});
        fprintf(fid,'\tPooled Data');
        for i=1:numel(pooledStatsNames)
            iStat = find(strcmp(pooledStatsNames{i},statsNames),1); 
            fprintf(fid,'\n%s\t',statsNames{iStat});
            fprintf(fid,'%g\t',statsData{iStat,:});
            fprintf(fid,'%g',pooledStatsData{i});
        end
        fclose(fid);
        
        if doPlot==1
            % save histograms of pooled distributions from iGroup
            plusTipMakeHistograms(groupData.M{iGroup},wtnDir);
            
            % make within-group boxplots (show each movie in iGroup)
            plusTipMakeBoxplots(groupData.dataMat{iGroup}',wtnGrpNames',wtnDir);
            
            % Plot comets as a function of time (Krek lab request)
            plotDetectedCometsNumber(groupData.detection{iGroup},wtnDir)
        end                  
    end
end

% Set up output directory
btwDir=[saveDir filesep 'btwGroupComparisons'];
if isdir(btwDir), rmdir(btwDir,'s'); end
mkdir(btwDir);

% Write stats results into a text file
statsFile = [btwDir filesep 'Stats.txt'];
statsNames = fieldnames(groupData.pooledStats{1});
pooledDataStats = cellfun(@struct2cell,groupData.pooledStats,'Unif',false);
pooledDataStats =horzcat(pooledDataStats{:});
fid=fopen(statsFile,'w+');
fprintf(fid,'\t%s',groupData.names{:});
for i=1:numel(statsNames)
    fprintf(fid,'\n%s\t',statsNames{i});
    fprintf(fid,'%g\t',pooledDataStats{i,:});
end
fclose(fid);

if doPlot
    % save histograms of pooled distributions from iGroup
    plusTipMakeHistograms(groupData.M,btwDir,'labels',groupData.names);
    
    % make between-group boxplots (show pooled data)
    pooledDataMat = cellfun(@(x) vertcat(x{:}),groupData.dataMat,'UniformOutput',false);
    plusTipMakeBoxplots(pooledDataMat',groupData.names',btwDir);
end


function plotDetectedCometsNumber(data,saveDir)

maxFrame=max(cellfun(@numel,data));
nProjects = numel(data);
colors=hsv(nProjects);

% define small and large fonts
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% plot
saveFig=figure('PaperPositionMode', 'auto'); % enable resizing
hold on;
arrayfun(@(i) plot(data{i},'-','Color',colors(i,:),'LineWidth',2),1:nProjects);

% Set thickness of axes, ticks and assign tick labels
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top','Xlim',[0 maxFrame]);
xlabel('Frame number', lfont{:});
ylabel('Comet number', lfont{:});

% remove box around figure
box off;

saveas(saveFig,[saveDir filesep 'DetectedCometsCount.tif']);
close(saveFig)