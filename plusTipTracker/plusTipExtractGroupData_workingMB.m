function [groupData]=plusTipExtractGroupData(groupList,varargin)
% plusTipPoolGroupData pools plus tip data from multiple projects in groups
%
% SYNOPSIS:  [groupData]=plusTipPoolGroupData(groupList,saveDir,doBtw,doWtn,doPlot,remBegEnd)
%
% INPUT:
% groupList : output of plusTipPickGroups, nProj x 2 cell array where the
%             first column contains the group identifier for each group,
%             and the second column contains the project path
% remBegEnd : 1 to remove tracks existing at the beginning
%             or end of the movie
%
% OUTPUT:
% groupData : structure containing group information and fields for 9
%             distributions: growth speed (gs), fgap speed (fs), bgap speed
%             (bs), growth lifetime (gl), fgap lifetime (fl), bgap lifetime
%             (bl), growth displacement (gd), fgap displacement (fd), and
%             bgap displacement (bd).
%
% Copyright (C) 2011 LCCB 
%
% This file is part of plusTipTracker.
% 
% plusTipTracker is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% plusTipTracker is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with plusTipTracker.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Maria Bagonis, April 2011
% Sebastien Besson, Seb 2011

%Input check
ip = inputParser;
ip.addRequired('groupList',@(x)iscell(x) || isempty(x));
ip.addOptional('remBegEnd',1,@isscalar);
ip.parse(groupList,varargin{:})
remBegEnd=ip.Results.remBegEnd;
if isempty(groupList), groupList=combineGroupListFiles; end


projGroupName=groupList(:,1);
projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'UniformOutput',0);

% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) ['grp_' regexprep(x,'[ -]','_')],...
    projGroupName,'uniformoutput',0);

% count unique groups and keep them in order of the original list
[btwGrpNames,m] = unique(projGroupName);
[~,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);
groupData.names = btwGrpNames';

M=cell(1,length(btwGrpNames));
S=cell(1,length(btwGrpNames));
Sgroup=cell(1,length(btwGrpNames));
Mgroup = cell(1,length(btwGrpNames)); 
D=cell(1,length(btwGrpNames));
dataByProject=cell(1,length(btwGrpNames));
for iGroup = 1:length(btwGrpNames)
    
    % indices of projects in iGroup
    projIndx=find(strcmp(btwGrpNames(iGroup),projGroupName));
    
    trkCount=1;
    for i = 1:length(projIndx)
        iProj = projIndx(i);
        
        % Read detection info
        s = load([projGroupDir{iProj} filesep 'feat' filesep 'movieInfo']);
        D{iGroup}{i,1}=arrayfun(@(x) size(x.xCoord,1),s.movieInfo);
        
        % Read post-tracking info 
        s = load([projGroupDir{iProj} filesep 'meta' filesep 'projData']);
        
     
        
        %
        if isfield(s.projData,'dataMatAllSubTracksConverted')
        dataMat = s.projData.mergedDataMatAllSubTracksConverted;
        else 
            dataMat = s.projData.dataMat_FOR_STATS;
        end 
        
        if remBegEnd==1
            dataMat = plusTipRemBegEnd(dataMat,s.projData); 
            % this output has data at beginning/end removed and units
            
            % already converted
            %[~,~,dataMat]=plusTipMergeSubtracks(s.projData);
        %else 
            % this output just gives merged tracks without converting units
            % or removing beginning/end data
            %dataMat=plusTipMergeSubtracks(s.projData);
            %dataMat(:,6)=dataMat(:,6).* s.projData.secPerFrame; % convert lifetimes to seconds
            %dataMat(:,7)=dataMat(:,7).*(s.projData.pixSizeNm/1000); % convert displacements to microns
        end
        
        dirByProj{iGroup}{i} = s.projData.anDir; 
        
        % reassign the track numbers so when combined from multiple projects they don't repeat
        trkIdx=unique(dataMat(:,1));
        dataMat(:,1)=swapMaskValues(dataMat(:,1),trkIdx,trkCount:trkCount+length(trkIdx)-1);
        trkCount=trkCount+length(trkIdx);
        
        % assign matrix to cell array
        dataByProject{iGroup}{i}=dataMat;
        [S{iGroup}{i},M{iGroup}{i}]= plusTipDynamParam(dataMat,s.projData,0,1);
    end
    [Sgroup{iGroup}]= plusTipDynamParam(vertcat(dataByProject{iGroup}{:}),s.projData,1,0);
    Mgroup{iGroup}= vertcat(M{iGroup}{:});  
    if isfield(Sgroup{iGroup},'nTracksSubRoi')
        
    % for now just recalculate pooled stats quick fix  
    Sgroup{iGroup}.stats.growth_speed_mean_INSIDE_REGION = nanmean(Mgroup{iGroup}(:,10)); 
    Sgroup{iGroup}.stats.growth_speed_median_INSIDE_REGION = nanmedian(Mgroup{iGroup}(:,10)); 
    Sgroup{iGroup}.stats.growth_lifetime_mean_INSIDE_REGION = nanmean(Mgroup{iGroup}(:,11)); 
    Sgroup{iGroup}.stats.growth_lifetime_median_INSIDE_REGION = nanmedian(Mgroup{iGroup}(:,11)); 
%     Sgroup{iGroup}.stats.polarCoordMeanOfAllSubtracks = nanmean(Mgroup{iGroup}(:,12)); 
%     Sgroup{iGroup}.stats.polarCoordMedianOfAllSubtracks = nanmedian(Mgroup{iGroup}(:,12)); 
%     Sgroup{iGroup}.stats.polarCoordMeanOfAllSubtracks_INSIDE_REGION = nanmean(Mgroup{iGroup}(:,13)); 
%     Sgroup{iGroup}.stats.polarCoordMedianOfAllSubtracks_INSIDE_REGION = nanmedian(Mgroup{iGroup}(:,13)); 
%   
    
    %also quickly sort by angle 
%     Sgroup{iGroup}.stats.perpAngleMeanGrowthSpeed = nanmean(Mgroup{iGroup}(:,14)); 
%     Sgroup{iGroup}.stats.parAngleMeanGrowthSpeed = nanmean(Mgroup{iGroup}(:,15)); 
%     
    
    
    end 
end
groupData.pooledStats = cellfun(@(x) x.stats,Sgroup,'UniformOutput',0);
groupData.pooledM = Mgroup;  
groupData.stats = cellfun(@(x) cellfun(@(y) y.stats,x,'Unif',0),S,'Unif',0);
groupData.dataMat=dataByProject;
groupData.M = M;
groupData.detection = D;
groupData.dirByProj = dirByProj; 
%save('groupData','groupData'); 