function [indStats,crossStats]=testTrajectories(trajectoryData,groupOrInd)
%
%
%
%

%============
% TEST INPUT
%============

% get info from input

% number of groups to compare
nGroups = length(trajectoryData(:));

% whether to compare individual stats
if nargin < 2 | isempty(groupOrInd)
    groupOrInd = 'group'
end
switch groupOrInd
    case 'group'
        if nGroups == 1
            error('to compare groups we need at least two sets of trajectories')
        end
        doIndividual = 0;
    case 'ind'
        % we find the number of individuals on the fly - later we can adapt
        % for preallocation etc.
        doIndividual = 1;
    otherwise
        error('wrong input for groupOrInd')
end

%==================

% do inds first, then group anyway, if there's more than one

% create colormap
logColormap = repeatEntries([1 1 1;...
                             0 0 0;...
                             0 0 1;...
                             0 1 1;...
                             0 1 0;...
                             1 1 0;...
                             1 0.5 0;...
                             1 0 0;...
                             0.5 0 0],...
                        [1,20,6,14,20,20,20,20,20]);
cLimits = [-0.05,7];

if doIndividual
    indStats(1:nGroups) = struct('compStruct',[],'probMat',[],'outlierTest',[],'nTrajectories',[]);
    crossStats(1:nGroups,1:nGroups) = struct('probMat',[]);
    for iGroup = 1:nGroups
        
        % loop through all the individual trajectories, collect all data
        % and then calculate ranksums
        
        % count trajectories
        nTrajectories = length(trajectoryData(iGroup).individualStatistics);
        % prepare list of data to compare
        compStruct(1:nTrajectories) = struct('growthSpeeds',[],'shrinkageSpeeds',[],...
            'growthTimes',[],'shrinkageTimes',[]);
        % prepare otuput
        [probMatGS,probMatSS,probMatGT,probMatST] = deal(ones(nTrajectories));
        
        
        for ti = 1:nTrajectories
            % read data. Calculate indices first, though
            [growthIdx,growthGroupIdx] = ...
                trajectoryAnalysisMainCalcStatsFindGroups(...
                trajectoryData(iGroup).individualStatistics(ti).dataListGroup,1);
            nGrowthGroups = size(growthGroupIdx,1);
            [shrinkageIdx,shrinkageGroupIdx] = ...
                trajectoryAnalysisMainCalcStatsFindGroups(...
                trajectoryData(iGroup).individualStatistics(ti).dataListGroup,2);
            nShrinkageGroups = size(shrinkageGroupIdx,1);
            
            % speeds. Speeds persisting over longer periods of time are
            % repeated per interval
            compStruct(ti).growthSpeeds = repeatEntries(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                growthIdx,4), diff(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                growthIdx,[1:2]),1,2))*60;
            compStruct(ti).shrinkageSpeeds = repeatEntries(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                shrinkageIdx,4), diff(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                shrinkageIdx,[1:2]),1,2))*60;
            
            % growth/ shrinkage group times
            
            compStruct(ti).growthTimes = diff(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                growthIdx,[1:2]),1,2);
%             compStruct(ti).growthTimes = zeros(nGrowthGroups,1);
%             for i=1:nGrowthGroups
%                 compStruct(ti).growthTimes(i) =...
%                     sum(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
%                     growthGroupIdx(i,1):growthGroupIdx(i,2),7));
%             end
            compStruct(ti).shrinkageTimes = diff(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
                shrinkageIdx,[1:2]),1,2);
%             compStruct(ti).shrinkageTimes = zeros(nShrinkageGroups,1);
%             for i=1:nShrinkageGroups
%                 compStruct(ti).shrinkageTimes(i) =...
%                     sum(trajectoryData(iGroup).individualStatistics(ti).dataListGroup(...
%                     shrinkageGroupIdx(i,1):shrinkageGroupIdx(i,2),7));
%             end
            
            for tj = ti-1:-1:1
                % compare
                [probMatGS(ti,tj),probMatGS(tj,ti)] = ...
                    deal(ranksum(compStruct(ti).growthSpeeds,compStruct(tj).growthSpeeds));
                [probMatSS(ti,tj),probMatSS(tj,ti)] = ...
                    deal(ranksum(compStruct(ti).shrinkageSpeeds,compStruct(tj).shrinkageSpeeds));
                [probMatGT(ti,tj),probMatGT(tj,ti)] = ...
                    deal(ranksum(compStruct(ti).growthTimes,compStruct(tj).growthTimes));
                [probMatST(ti,tj),probMatST(tj,ti)] = ...
                    deal(ranksum(compStruct(ti).shrinkageTimes,compStruct(tj).shrinkageTimes));
                    
            end
        end
        
        % test for outliers
        outlierTest = zeros(nTrajectories,4);
        for ti = 1:nTrajectories
            % find the values of all except the one we want to compare
            allOtherGS = cat(1,compStruct([1:ti-1,ti+1:end]).growthSpeeds);
            allOtherSS = cat(1,compStruct([1:ti-1,ti+1:end]).shrinkageSpeeds);
            allOtherGT = cat(1,compStruct([1:ti-1,ti+1:end]).growthTimes);
            allOtherST = cat(1,compStruct([1:ti-1,ti+1:end]).shrinkageTimes);
            
            % compare
            outlierTest(ti,1)=ranksum(allOtherGS,compStruct(ti).growthSpeeds);
            outlierTest(ti,2)=ranksum(allOtherSS,compStruct(ti).shrinkageSpeeds);
            outlierTest(ti,3)=ranksum(allOtherGT,compStruct(ti).growthTimes);
            outlierTest(ti,4)=ranksum(allOtherST,compStruct(ti).shrinkageTimes);
            
        end
            
        
        % display
        dispMat1 = [probMatGS, repmat(1.122,nTrajectories,1), probMatSS; ...
                repmat(1.122,1,2*nTrajectories+1);...
                probMatGT, repmat(1.122,nTrajectories,1),probMatST];
        dispMatLog = -log10(dispMat1);
        uH=uiviewpanel;
        imshow(dispMatLog);
        % apply colormap and cLims
        set(uH,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        
        
        
        
        % outliers
        uiviewpanel;
        outlierTestLog = -log10(outlierTest);
        imshow(outlierTestLog);
        % apply colormap and cLims
        set(gcf,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        
        
        % store data
        indStats(iGroup).compStruct = compStruct;
        indStats(iGroup).probMat = cat(3,probMatGS, probMatSS, probMatGT, probMatST);
        indStats(iGroup).outlierTest = outlierTest;
        indStats(iGroup).nTrajectories = nTrajectories;
        
        % cross-test between groups
        for jGroup = iGroup-1:-1:1
            [probMatGS,probMatSS,probMatGT,probMatST] = ...
                deal(ones(indStats(iGroup).nTrajectories,indStats(jGroup).nTrajectories));
            for ti = 1:indStats(iGroup).nTrajectories
                % loop through full group J, because we won't get anything
                % symmetric
                for tj = 1:indStats(jGroup).nTrajectories
                    probMatGS(ti,tj)=ranksum(indStats(iGroup).compStruct(ti).growthSpeeds,indStats(jGroup).compStruct(tj).growthSpeeds);
                    probMatSS(ti,tj)=ranksum(indStats(iGroup).compStruct(ti).shrinkageSpeeds,indStats(jGroup).compStruct(tj).shrinkageSpeeds);
                    probMatGT(ti,tj)=ranksum(indStats(iGroup).compStruct(ti).growthTimes,indStats(jGroup).compStruct(tj).growthTimes);
                    probMatST(ti,tj)=ranksum(indStats(iGroup).compStruct(ti).shrinkageTimes,indStats(jGroup).compStruct(tj).shrinkageTimes);
                end
            end
               
            crossStats(iGroup,jGroup).probMat = cat(3,probMatGS, probMatSS, probMatGT, probMatST);
        end
        
        
    end %for iGroup = 1:nGroups
    
    % display big matrices
    nTrajectoriesAll = cat(1,indStats.nTrajectories);
    % allocate matrices. We want separations between the individual
    % matrices. Remember that we take the log of the thing
    bigMatSize = sum(nTrajectoriesAll) + length(nTrajectoriesAll)-1;
    [allGS,allSS,allGT,allST] = deal(repmat(1.122,bigMatSize,bigMatSize));
    currentRow = 0;
   
    
    for iGroup = 1:nGroups % rows
        currentCol = 0;
        for jGroup = 1:nGroups % cols - start at zero
            % decide from where we want to read the information
            switch (iGroup == jGroup) + 2*(iGroup < jGroup)
                % 0: below diagonal
                % 1: diagonal
                % 2: above diagonal
                case 1
                    % read from indStats
                    allGS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        indStats(iGroup).probMat(:,:,1);
                    allSS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        indStats(iGroup).probMat(:,:,2);
                    allGT(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        indStats(iGroup).probMat(:,:,3);
                    allST(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        indStats(iGroup).probMat(:,:,4);
                    
                case 0
                    % read from crossStats, don't transpose
                    allGS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(iGroup,jGroup).probMat(:,:,1);
                    allSS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(iGroup,jGroup).probMat(:,:,2);
                    allGT(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(iGroup,jGroup).probMat(:,:,3);
                    allST(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(iGroup,jGroup).probMat(:,:,4);
                    
                case 2
                    % read from crossStats, transpose
                    allGS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(jGroup,iGroup).probMat(:,:,1)';
                    allSS(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(jGroup,iGroup).probMat(:,:,2)';
                    allGT(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(jGroup,iGroup).probMat(:,:,3)';
                    allST(currentRow+1:currentRow+nTrajectoriesAll(iGroup),currentCol+1:currentCol+nTrajectoriesAll(jGroup))=...
                        crossStats(jGroup,iGroup).probMat(:,:,4)';
            end % switch
            currentCol = currentCol+nTrajectoriesAll(jGroup)+1;
        end % for jGroup = 1:nGroups - cols
        currentRow = currentRow+nTrajectoriesAll(iGroup)+1;
    end % for iGroup = 1:nGroups - rows
    
    % display
    uH=figure('Name','GrowthSpeed');
        imshow(-log10(allGS));
        % apply colormap and cLims
        set(uH,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        colorbar
        
        uH=figure('Name','ShrinkageSpeed');
        imshow(-log10(allSS));
        % apply colormap and cLims
        set(uH,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        colorbar
        
        uH=figure('Name','Catastrophe');
        imshow(-log10(allGT));
        % apply colormap and cLims
        set(uH,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        colorbar
        
        uH=figure('Name','Rescue');
        imshow(-log10(allST));
        % apply colormap and cLims
        set(uH,'Colormap',logColormap);
        set(gca,'CLim',cLimits)
        colorbar
    
end %if doIndividual

%====================
% COMPARE GLOBALLY
%====================

% read values
if doIndividual
    % read global params from indStats
    for i = 1:nGroups
        cs = indStats(i).compStruct;
        globStats(i).growthSpeeds = catStruct(1,'cs.growthSpeeds');
        globStats(i).shrinkageSpeeds = catStruct(1,'cs.shrinkageSpeeds');
        globStats(i).growthTimes = catStruct(1,'cs.growthTimes');
        globStats(i).shrinkageTimes = catStruct(1,'cs.shrinkageTimes');
    end
else
end

% compare all
[globCompareGS,globCompareSS,globCompareGT,globCompareST] = deal(zeros(nGroups));
for gi = 1:nGroups
    for gj = gi:-1:1
        [globCompareGS(gi,gj),globCompareGS(gj,gi)] = deal(...
            ranksum(globStats(gi).growthSpeeds,globStats(gj).growthSpeeds));
        [globCompareSS(gi,gj),globCompareSS(gj,gi)] = deal(...
            ranksum(globStats(gi).shrinkageSpeeds,globStats(gj).shrinkageSpeeds));
        [globCompareGT(gi,gj),globCompareGT(gj,gi)] = deal(...
            ranksum(globStats(gi).growthTimes,globStats(gj).growthTimes));
        [globCompareST(gi,gj),globCompareST(gj,gi)] = deal(...
            ranksum(globStats(gi).shrinkageTimes,globStats(gj).shrinkageTimes));
    end
end

% display
uH = uiviewpanel;
set(uH,'Name','GlobalGrowth')
imshow(-log10(globCompareGS));
set(uH,'Colormap',logColormap);
set(gca,'CLim',cLimits)
uH = uiviewpanel;
set(uH,'Name','GlobalShrinkage')
imshow(-log10(globCompareSS));
set(uH,'Colormap',logColormap);
set(gca,'CLim',cLimits)
uH = uiviewpanel;
set(uH,'Name','GlobalCatastrophe')
imshow(-log10(globCompareGT));
set(uH,'Colormap',logColormap);
set(gca,'CLim',cLimits)
uH = uiviewpanel;
set(uH,'Name','GlobalRescue')
imshow(-log10(globCompareGS));
set(uH,'Colormap',logColormap);
set(gca,'CLim',cLimits)
