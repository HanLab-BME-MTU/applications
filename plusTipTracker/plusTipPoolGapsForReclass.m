function [projData ] = plusTipPoolGapsForReclass( groupList, varargin)
%If Unimodal thresholding requires a full population of fgaps 
%  (ie the number of fgaps per cell is too small to get reliable
%  thresholding. One can use this function to pool fgap data and 
% apply the threshold. 
%% 
% Input 
% groupList:  groupList of projects of the same condition
%
% meta2Use: meta folder from which you would like to pool your data
% (typcically  'meta') input as a string
%
% saveCopy: 1 if you want to save a copy of the old projData before the
% reclass
%
% metaOldNewName: if saveCopy = 1, give a name to the folder (if saveCopy =0 
% just set to []) 
%
% mkHistGaps: set to 1 to save a copy of the histogram of the pooled fgap 
% speeds including the unimodal thresh  
% 
% useFirstInList: set to 1 to use the first group in list to calculate the 
% threshold for all groups in the groupList
% 
% mkHist: set to 1 to save histograms of speed, lifetime, displacement 
% (same option as in the post-processing settings)

ip = inputParser;
ip.addRequired('groupList',@(x)iscell(x) || isempty(x));
if nargin<1 || isempty(groupList)
    [file2load dir] = uigetfile(pwd, 'Please Choose a groupList file to Load');
    s = load([dir filesep file2load]); 
    groupList= s.groupList;
end 

ip.addOptional('groupListForThresh',groupList,@iscell);
ip.addOptional('meta2Use','meta',@ischar);
ip.addOptional('saveCopy',1,@isscalar);
ip.addOptional('oldMetaSaveName','metaBeforeReclass',@ischar);


%ip.addParam('groupListForThresh',groupList, @(x)iscell(x) || isempty(x));
ip.addParamValue('mkHistGaps',1,@isscalar);
ip.addParamValue('useFirstInList',1,@isscalar);
ip.addParamValue('mkHist',1,@isscalar);
ip.addParamValue('remBegEnd',1,@isscalar);


ip.parse(groupList,varargin{:});

meta2Use = ip.Results.meta2Use;
saveCopy = ip.Results.saveCopy;
oldMetaSaveName = ip.Results.oldMetaSaveName;
groupListForThresh = ip.Results.groupListForThresh;
mkHistGaps = ip.Results.mkHistGaps;
useFirstInList = ip.Results.useFirstInList;
mkHist = ip.Results.mkHist;
remBegEnd = ip.Results.remBegEnd; 



%% Choices for Input  (FIX)
onlyFgap2GrowthReclass = 1; %  unimodal thresholding
% on pooled gap data


bgapUniModeThreshNoCorrect = 0; % use symmetrical reclassification (currenlty default is to not reclassify bgaps as pauses)

bgapUniModeThreshCorrect = 0; % attempts to correct for the overestimation of the bgap to pause threshold that occurs due to comet latency

bgapReclassFluctRadius = 0; %   not recommended



%% collect fgaps for thresholding


projGroupNameThresh=groupListForThresh(:,1);
projGroupDirThresh= groupListForThresh(:,2);



% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupNameThresh=cellfun(@(x) strrep(x,'-','_'),projGroupNameThresh,'uniformoutput',0);
projGroupNameThresh=cellfun(@(x) strrep(x,' ','_'),projGroupNameThresh,'uniformoutput',0);
projGroupNameThresh=cellfun(@(x) ['grp_' x],projGroupNameThresh,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupNameThresh);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);

groupNames = cell(length(btwGrpNames),1);
for iGroup = 1:length(btwGrpNames)
    groupNames(iGroup,1) = btwGrpNames(iGroup,1);
end


for iGroup = 1:length(btwGrpNames)
    
    
    
    tempIdxThresh= strmatch(btwGrpNames(iGroup),projGroupNameThresh,'exact');
    
    %if (iGroup > 1 && useFirstInList == 1); % do nothing already calculated cut-offs
        
        
            
        
   % else % recalculate cut-off based on each group
        
   % collect all data
        for iProj = 1:length(tempIdxThresh)
            dirToload = formatPath([projGroupDirThresh{tempIdxThresh(iProj)} filesep meta2Use ]);
            temp = load([dirToload filesep 'projData.mat']);
            
            
            projData = temp.projData;
            
            dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
            
            if iProj == 1
                dataMatPooled = dataMat;
            else
                dataMatPooled = [dataMatPooled ; dataMat] ;
            end
            
        end % for iProj
        
         
        
        % perform unimodal thresholding using all fgap speeds to obtain
        % maximum fgap speed corresponding to a pause event (above this thresh
        % are likely undetected growth events)
        
        fgapIdx = find(dataMatPooled(:,5) == 2 | dataMatPooled(:,5) == 5);
        
        
        
        fgapSpeeds = dataMatPooled(fgapIdx,4);
        
        %gapSpeed95th = prctile(fgapSpeeds,95);
        
        
        %fgapSpeeds = fgapSpeeds(fgapSpeeds<fgapSpeed95th);
        
        if (iGroup > 1 && useFirstInList == 1) % skip recalculation of cut-off
        else % recalculate
            
            [cutoffIdx, cutoffValueFGap, sp, axesH,maxBinValue] = cutFirstHistMode(fgapSpeeds,mkHistGaps);
            
          
            cutOffValueFGap_VelMicPerMin = cutoffValueFGap;
            maxBinValue_VelMicPerMin = maxBinValue;
            numFgapsTotal = length(fgapIdx);
            fgapReclassScheme = 'Unimodal Thresholding: Pooled Data';
            if mkHistGaps == 1
                [pathup1] =  getFilenameBody(projGroupDirThresh{tempIdxThresh(iProj)});
                [pathup2 name ] = getFilenameBody(pathup1);
                pathup2 = formatPath(pathup2);
                
                uniModalFigureDir = [pathup2 filesep 'PoolUnimodalThreshFigures'];
                
                if ~isdir(uniModalFigureDir)
                    mkdir(uniModalFigureDir);
                end
                saveas(gcf,[uniModalFigureDir filesep name 'UnimodalThresh.fig']);
                close(gcf);
            end
        end % iGroup>1 && useFirstInList
        
        
        % always make gap speed histograms
        if mkHistGaps == 1
            
           
            
           
      
            bgapIdx = dataMatPooled(:,5) ==3 | dataMatPooled(:,5)==6; 
            bgapSpeeds = dataMatPooled(bgapIdx,4);  
            allGapSpeeds = [fgapSpeeds;bgapSpeeds] ; 
            
            binCenters = floor(min(allGapSpeeds)):1:ceil(max(allGapSpeeds)); 
            h = figure;
            [n, xout] = hist(allGapSpeeds,binCenters); 
            bar(binCenters,n/sum(n)); 
            hold on 
            % line to mark threshold 
            x = ones(1,max(n)+10+1); 
            x = x*cutOffValueFGap_VelMicPerMin; 
            y= (0:max(n)+10); 
            xzeros = zeros(1,max(n)+10+1);
            plot(x,y,'r'); 
            plot(xzeros,y,'r'); 
            
            forPlotNames = cellfun(@(x) strrep(x,'_',''),btwGrpNames(:,1),'uniformoutput',0); 
            title([forPlotNames(iGroup,1),'Velocity All Gaps']); 
            xlabel('Gap Velocity (um*min-1)');
            ylabel('Number Of Gaps');
            axis([min(allGapSpeeds) max(allGapSpeeds) 0 0.1]);%(max(n)+10)
            filename = [char(btwGrpNames(iGroup,1)), 'gapHist'];
            saveas(gcf,[uniModalFigureDir filesep filename '.eps'],'psc2');
            if ~isdir([uniModalFigureDir filesep 'tifs']) 
                mkdir([uniModalFigureDir filesep 'tifs']);
            end 
            saveas(gcf,[uniModalFigureDir filesep 'tifs' filesep filename '.tif']); 
            saveas(gcf,[uniModalFigureDir filesep filename,'fig'],'fig'); 
            
        end % if mkGapHist
        %% BGap2Pause No Reclass
        if onlyFgap2GrowthReclass == 1 %
            bgapThresh = 0;
            bgapReclassScheme = 'No Bgap2Pause Reclass';
            cutOffValueBGap_VelMicPerMin = 0;
        end
        %%  Bgap Reclassification Using Unimodal Fgap Thresh
        
        if bgapUniModeThreshNoCorrect == 1
            
            bgapThresh = cutoffValueFGap;
            bgapReclassScheme = 'Unimodal Fgap Thresh- No Correct: Pooled Data';
            cutOffValueBGap_VelMicPerMin = cutoffValueFGap;
            
        end
        
        %%  Bgap Reclassification Using "Corrected" Unimodal Fgap Thesh
        
        if bgapUniModeThreshCorrect == 1
            
            widthDist = cutoffValueFGap - maxBinValue;
            
            
            
            if widthDist < maxBinValue % distribution does not go into the negative values
                bgapThresh = 0;
            else
                bgapThresh = widthDist - maxBinValue;
            end
            
            bgapReclassScheme = 'Unimodal Fgap Thresh- Correct For Comet Latency: Pooled Data';
            cutOffValueBGap_VelMicPerMin = bgapThresh;
            maxBinValue_VelMicPerMin = maxBinValue;
            fgapDisbWidth_VelMicPerMin = widthDist;
            
            
        end
        
    
    
    %% Using Above Thresholds Correct Individual Projects
    
    % set up directories if different from
    projGroupName=groupList(:,1);
    projGroupDir= groupList(:,2);
    
    
    
    % fix the names if there are spaces or hyphens and append prefix 'grp'
    projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
    projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
    projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);
    
    
    % count unique groups and keep them in order of the original list
    [btwGrpNames,m,projGroupIdx] = unique(projGroupName);
    [b,idx]=sort(m);
    btwGrpNames=btwGrpNames(idx);
    
     tempIdx= strmatch(btwGrpNames(iGroup),projGroupName,'exact');
    
    for iProj = 1:length(tempIdx)
        dir2Load = formatPath([projGroupDir{tempIdx(iProj)} filesep meta2Use]); % potentially use separate dir
        temp = load([dir2Load filesep 'projData']);
        projData = temp.projData;
        dataMat =  projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
        
        
        fgapIdx = find(dataMat(:,5) == 2 | dataMat(:,5) == 5);
        bgapAllIdx = find(dataMat(:,5) == 3 | dataMat(:,5) == 6);
        
        
        dataMat(fgapIdx,5) = 2;  % make sure to undo any previous reclass scheme
        dataMat(bgapAllIdx,5) = 3;
        
        
        growthFgapIdx = find((dataMat(:,5) == 2 | dataMat(:,5) == 5) & abs(dataMat(:,4)) > cutoffValueFGap);
        dataMatReclass = dataMat;
        dataMatReclass(growthFgapIdx,5) =  5;
        
        if bgapReclassFluctRadius == 1
            bgap2pauseIdx = find((dataMat(:,5) == 3 & abs(dataMat(:,7)) < projData.trackingParameters.fluctRadius));
        else
            
            
            bgap2pauseIdx = find((dataMat(:,5) == 3 | dataMat(:,5) == 6) & abs(dataMat(:,4)) < bgapThresh);
            dataMatReclass(bgap2pauseIdx,5) = 6;
            dataMat(bgap2pauseIdx,5) = 2;
            projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix = dataMatReclass;
        end
        
        %% SAVE INFO REGARDING THE RECLASS SCHEME EMPLOYED %%
        if useFirstInList == 1
            
            projData.fgapReclassScheme = [fgapReclassScheme ': Use ' char(btwGrpNames(1)) ' To Set Thresh '];
            
            projData.bgapReclassScheme = [bgapReclassScheme ': Use ' char(btwGrpNames(1)) ' To Set Thresh '];
            
        else
            projData.fgapReclassScheme = fgapReclassScheme;
            projData.bgapReclassScheme = bgapReclassScheme;
        end
        
        projData.cutOffValueFGap_VelMicPerMin = cutOffValueFGap_VelMicPerMin;
        projData.numFGapsUniModalHist = numFgapsTotal;
        
        if (bgapUniModeThreshCorrect == 1 || bgapUniModeThreshNoCorrect ==1 ||onlyFgap2GrowthReclass == 1 )
            projData.cutOffValueBGap_VelMicPerMin= cutOffValueBGap_VelMicPerMin;
        end
        
        if bgapUniModeThreshCorrect == 1
            projData.maxBinValue_VelMicPerMin = maxBinValue_VelMicPerMin;
            projData.fgapDisbWidth_VelMicPerMin = fgapDisbWidth_VelMicPerMin;
        end
        
        
        projData.tracksWithFGapPause = dataMatReclass(dataMatReclass(:,5) == 2,1);
        projData.tracksWithFGapGrowth = dataMatReclass(growthFgapIdx,1);
        projData.tracksWithBGapPause = dataMatReclass(bgap2pauseIdx,1);
        projData.tracksWithBGapShrink = dataMatReclass(dataMatReclass(:,5) == 3,1);
        
        
        
        %%
        
        
        % MERGE RECLASSIFICATIONS (ie subTrack ID 5 --> 1 need to incorporate reclassified pauses in growth stats)
        % merging will be performed regardless of reclassification scheme unless we
        % are pooling data
        
        
        % these are the affected track numbers for growth fgaps
        tracks2check=unique(dataMat(growthFgapIdx,1));
        rows2remove=[];
        
        for i=1:length(tracks2check)
            % rows of dataMat corresponding to i track
            subIdx=find(dataMat(:,1)==tracks2check(i));
            % rows of dataMat corresponding to fgaps in track to consolidate
            fgap2remIdx=intersect(growthFgapIdx,subIdx);
            % rows of dataMat corresonding to bgaps or real pauses in track
            sepIdx=union(subIdx(dataMat(subIdx,5)==3),setdiff(intersect(subIdx,fgapIdx),fgap2remIdx))';
            % split the track based on bgaps, so that all fgaps that
            % should be consolidated together can be done dataMat the same time
            sIdx=[subIdx(1); sepIdx(:)];
            eIdx=[sepIdx(:); subIdx(end)];
            % loop through groups of subtracks (split by bgaps)
            for j=1:length(sIdx)
                % pTemp contains the ones to consolidate in this section
                pTemp=intersect(fgap2remIdx,sIdx(j):eIdx(j));
                if ~isempty(pTemp)
                    fIdx=min(pTemp)-1; % first row - prior to first fgap
                    lIdx=max(pTemp)+1; % last row - after final fgap
                    
                    dataMat(fIdx,3)=dataMat(lIdx,3); % change end frame
                    dataMat(fIdx,6)=sum(dataMat(fIdx:lIdx,6)); % sum lifetimes
                    dataMat(fIdx,7)=sum(dataMat(fIdx:lIdx,7)); % sum total displacements
                    dataMat(fIdx,4)= pixPerFrame2umPerMin(dataMat(fIdx,7)/dataMat(fIdx,6),projData.secPerFrame,projData.pixSizeNm); % find new average velocity
                    
                    % keep track of which are the extra rows
                    rows2remove=[rows2remove fIdx+1:lIdx];
                end
            end
        end
        
        % remove the extra rows
        dataMat(rows2remove,:)=[];
        
        % Conversions and Calculation of Comet Latency
        % calculate necessary avg values from above merged data (where all
        % reclassified pauses have been merged with preceding growth
        % subtrack
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add 2 more columns corresponding to extrac subTrack information
        % Though obvious here when partition these subtracks based on
        % regional criteria it can become ambigious.
        % Therefore it is useful to mark these here while the dataStruct for the
        % whole cell compound tracks are in tact.
        % Column 8: Nuc Event 1= yes 0 = no
        % Column 9: Growth Before Term Event = 1, Growth Before Fgap = 2,
        % Growth before Bgap = 3, Growth before undefined gap = 4
        
        [ dummy nucEventsIdx] =  unique(dataMat(:,1),'first');
        [ dummy termEventsIdx dummy] = unique(dataMat(:,1),'last');
        
        dataMat(nucEventsIdx,8) = 1;
        allIdx = 1:length(dataMat(:,1));
        nonNucIdx = setdiff(allIdx,nucEventsIdx);
        
        dataMat(nonNucIdx,8) = 0;
        
        dataMat(termEventsIdx,9) = 1;
        
        fIdx = find(dataMat(:,5) == 2);
        bIdx = find(dataMat(:,5) == 3);
        uIdx = find(dataMat(:,5) == 4);
        
        beforeFgapIdx=fIdx-1;
        beforeBgapIdx = bIdx-1;
        beforeUgapIdx = uIdx-1;
        
        dataMat(beforeFgapIdx,9) = 2;
        dataMat(beforeBgapIdx,9) = 3;
        dataMat(beforeUgapIdx,9) = 4;
        dataMat([fIdx;bIdx;uIdx],9)= 0;
        
        
        
        % perform lifetime and displacement unit conversions
        dataMat(:,6)=dataMat(:,6).* projData.secPerFrame; % convert lifetimes to seconds
        dataMat(:,7)=dataMat(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns
        
        % mark if part of compound versus non-compound track
        gapIdx = sort([fIdx;bIdx]);
        compIdx= unique(sort([gapIdx ; (gapIdx+1) ; (gapIdx -1)]));
        % set those part of a compound track to in column 10 to 1
        % so marked for later partitioning
        dataMat(compIdx,10) = 1;
        
        compDataMat = dataMat(compIdx,:);
        
        
        
        % Segregate Tracks That Are Exclusively From Single Tracks
        singleDataMat = dataMat;
        singleDataMat(compIdx,:) = [];
        
        % remove uIdx
        uIdxSingleMat = find(singleDataMat(:,5) == 4);
        toRemove = sort([uIdxSingleMat;uIdxSingleMat+1;uIdxSingleMat-1]);
        toRemove = toRemove(toRemove~=0);
        singleDataMat(toRemove,:) = [];
        
        
        if remBegEnd == 1
            compDataMat = plusTipRemBegEnd(compDataMat,projData,1);
            singleDataMat = plusTipRemBegEnd(singleDataMat,projData,1);
        end
        
        projData.compDataMat = compDataMat;
        projData.singleDataMat = singleDataMat;
        
        % save the merged tracks dataStructure in projData
        projData.mergedDataMatAllSubTracksConverted = dataMat; % save this dataStruct for subRoi
        % partitioning and stat calculations from pooled data.
        
        
        
        % Remove growth subtracks that start in the first frame and end in the
        % last frame (as well as flanking fgap and bgaps)
        if remBegEnd == 1
            dataMatCrpSecMic = plusTipRemBegEnd(dataMat,projData);
            projData.remBegEnd = 'yes';
        else
            projData.remBegEnd = 'no';
        end
        
        %% Calculate Stats from the matrix where the tracks at the beginning and end have been removed
        [projData,projData.M] = plusTipDynamParam(dataMatCrpSecMic,projData,0,0);
        
        MCell = num2cell(projData.M);
        titles = cell(1,9);
        titles{1,1} = 'growth speed (microns/min)';
        titles{1,2} = 'fgap speed (microns/min)';
        titles{1,3} = 'bgap speed (microns/min)';
        titles{1,4} = 'growth lifetimes (sec)';
        titles{1,5} = 'fgap lifetimes (sec)';
        titles{1,6} = 'bgap lifetimes (sec)';
        titles{1,7} = 'growth displacements (microns)';
        titles{1,8} = 'fgap displacements (microns)';
        titles{1,9} = 'bgap displacementes (microns)';
        
        
        projData.MCell = [titles;MCell];
        
        %avgVelGrowth = mean(dataMat(dataMat(:,5) == 1,4));
        
        %avgDispPause = mean(dataMat(dataMat(:,5) == 2,7)); % in microns
        %absAvgDispPause = mean(abs(dataMat(dataMat(:,5) == 2,7))); % in microns consider bgaps positive
        %avgDispPauseNoBgapReclass= mean(dataMatReclass(dataMatReclass(:,5) == 2, 7)).*projData.pixSizeNm/1000;
        
        
        % calc avg latency of comet formation  (in min)
        %avgLat = avgDispPause/avgVelGrowth; % in minutes
        %projData.avgCometLatSec = avgLat*60; % in seconds
        
        %absAvgLat = absAvgDispPause/avgVelGrowth;
        %projData.avgCometLatSecAbs = absAvgLat*60;
        
        %avgLatNoBgapReclass = avgDispPauseNoBgapReclass/avgVelGrowth; % in minutes
        %projData.avgCometLatSecNoBgapReclass = avgLatNoBgapReclass*60;
        
        
        %Calculate Fraction of SubTracks Reclassified
        % fraction of bgaps that are changed to pause, because their speeds are
        % slower than the cutoff for fgap pauses
        
        if isempty(bgapAllIdx)
            projData.percentBgapsReclass=NaN;
        else
            projData.percentBgapsReclass=100*length(bgap2pauseIdx)/length(bgapAllIdx);
        end
        
        % fraction of fgaps that were consolidated into growth, since their speeds
        % were more than 50% of the growth speed just prior to the gap
        if isempty(fgapIdx)
            projData.percentFgapsReclass=NaN;
        else
            projData.percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
        end
        if saveCopy  == 1
            % save a copy of the meta Dir Before Pooling
            dirBeforeRename = dir2Load;
            dirAfterRename = regexprep(dir2Load,meta2Use,oldMetaSaveName);
            
            if isdir(dirAfterRename)
                rmdir(dirAfterRename,'s');
            end
            
            movefile(dirBeforeRename, dirAfterRename);
        else % don't save a copy
            % rmdir(dirBeforeRename)
        end
        
        
        % rewrite the projData
        dir2Write = regexprep(dir2Load,meta2Use,'meta');
        if isdir(dir2Write) == 0
            mkdir(dir2Write);
        end
        save([dir2Write filesep 'projData'],'projData')
        
        
        if mkHist==1
            plusTipMakeHistograms(projData.M,[dir2Write filesep 'histograms'])
            %plusTipPlotTrackAngles(projData,[dir2Write filesep 'histograms']);
        end
        
        
    end % for iProj
end % for iGroup
end
