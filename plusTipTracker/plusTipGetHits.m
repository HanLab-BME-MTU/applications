function [dataStruct] = plusTipGetHits(saveDir,groupList,runStats,getHits,getMedValues,stringency,testID1, testID2)


% Extracts Different Stats from meta data stored in projData.stats (Both
% cell autonomous parameters (ie % growth terminal etc). AND mean
% and med perTrackParameters (ie growth speed, etc) for each cell.  
% NOTE: Assumes each project is an individual cell 
% Groups these values for each cell (ie each project) based on groupList and feeds into the discrimMat 
% (can perform either the permutation T-test or t-test to validate if the two populations 
% are equivalent. 
% can read out hits 
% NOTE: Certainly NOT optimally written: Quick Fix MB 03/11
% 
% 
% Input: 
% saveDir = directory to which you would like to save the hits, if [] it
% will ask you to specify a directory
% 
% groupList = output of createGroupList in Gui (just load by double
% clicking)
%
% if runStats = 1: run the discrimination matrix
% if getHits = 1: output hits from discrimination matrix based on cut-off
%                 (currently reads hits from the first column of
%                 discrimMat)
%
% stringency = P-value cut-off for defining hits
% 
% testID1 = test below diagonal in discrimMat: test for which hits will be found  
%
% testID2 = test above diagonal in discrimMat: statistical values for this test will be 
%           recorded and save in the discrimMat but hits will NOT be read out. 
% IDs for tests that can be read performed by disciminationMatrix.m
% 
%                          1 - t-test for means (default below diag.)
%                          2 - Wilcoxon ranksum test for medians
%                          10 - K-S test for distributions
%                          11 - K-S test for distributions with subtracted
%                               mean (default for above diag.)
%                          12 - K-S test for distributions with subtracted
%                               median
%                          20 - permutation test for means
%                          21 - distribution test - calibrated K-S test
%                               with mean subtraction
%


if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for hit files.');
end

    
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

groupNames = cell(length(btwGrpNames),1);
for iGroup = 1:length(btwGrpNames)
    groupNames(iGroup,1) = btwGrpNames(iGroup,1);
end 


for iGroup = 1:length(btwGrpNames)

    % indices of projects in iGroup
    tempIdx=strmatch(btwGrpNames(iGroup),projGroupName,'exact');
   
    %reset
    
    gFor = zeros(length(tempIdx),1);
    gBack = zeros(length(tempIdx),1);
    gTerm = zeros(length(tempIdx),1); 
    ratioPause2Growth = zeros(length(tempIdx),1);
    ratioShrink2Growth = zeros(length(tempIdx),1);
    dynamicity = zeros(length(tempIdx),1);
    ratioTotalTrack2GrowthSubtrack = zeros(length(tempIdx),1);
    ratioGapTracks2NonGapTracks = zeros(length(tempIdx),1);
    cometDensity = zeros(length(tempIdx),1); 
    avgCometLat = zeros(length(tempIdx),1);
    avgCometLatBeforeReclass = zeros(length(tempIdx),1);
    fgap_freqTime = zeros(length(tempIdx),1);
    fgap_freqLength = zeros(length(tempIdx),1);
    bgap_freqTime = zeros(length(tempIdx),1);
    bgap_freqLength = zeros(length(tempIdx),1);
    term_freqTime = zeros(length(tempIdx),1);
    term_freqLength = zeros(length(tempIdx),1);
    freqTimeFgap2freqTimeTerm  = zeros(length(tempIdx),1);
    freqLengthFgap2freqLengthTerm = zeros(length(tempIdx),1);
    freqTimeBgap2freqTimeFgap = zeros(length(tempIdx),1);
    
    gsAvgPerCell = zeros(length(tempIdx),1);
    glAvgPerCell = zeros(length(tempIdx),1);
    gdAvgPerCell = zeros(length(tempIdx),1);
    fsAvgPerCell = zeros(length(tempIdx),1);
    flAvgPerCell = zeros(length(tempIdx),1);
    fdAvgPerCell = zeros(length(tempIdx),1);
    bsAvgPerCell = zeros(length(tempIdx),1);
    blAvgPerCell = zeros(length(tempIdx),1);
    bdAvgPerCell = zeros(length(tempIdx),1);
    
    
    
    
    if getMedValues == 1
    %Median Values:
 
    gsMedPerCell = zeros(length(tempIdx),1);
    glMedPerCell = zeros(length(tempIdx),1);
    gdMedPerCell = zeros(length(tempIdx),1);
    fsMedPerCell = zeros(length(tempIdx),1);
    flMedPerCell = zeros(length(tempIdx),1);
    fdMedPerCell = zeros(length(tempIdx),1);
    bsMedPerCell = zeros(length(tempIdx),1);
    blMedPerCell = zeros(length(tempIdx),1);
    bdMedPerCell = zeros(length(tempIdx),1);
    else 
    end % getMedValues
    
    dataByProject=cell(length(tempIdx),1);
    trkCount=1;
    
    
    % Collect the Avg Value From Each Project in the Group
    
    for iProj = 1:length(tempIdx)
        temp = load([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData']);
        gFor(iProj,1) = temp.projData.stats.percentGrowthLinkedForward;
        gBack(iProj,1) = temp.projData.stats.percentGrowthLinkedBackward;
        gTerm(iProj,1) = temp.projData.stats.percentGrowthTerminal;
        ratioPause2Growth(iProj,1) = temp.projData.stats.nFgaps/temp.projData.stats.nGrowths;
        ratioShrink2Growth(iProj,1) = temp.projData.stats.nBgaps/temp.projData.stats.nGrowths;
        dynamicity(iProj,1) = temp.projData.stats.dynamicity;
        
        trackInfo = temp.projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
        
        numGrowthSubTracksAll = length(trackInfo(trackInfo == 1));
        ratioTotalTrack2GrowthSubtrack(iProj,1) = temp.projData.nTracks/numGrowthSubTracksAll; 
        %  (pooled over a cell)
        % values closer to 1 suggest less fragmentation. note this parameter is calculated for
        % the ENTIRE POPULATION Where tracks initiated in the first and
        % last frame of the movie included
         
        
        idxTracksWithGap = find(trackInfo(trackInfo(:,5) == 2 ...
        | trackInfo(:,5) == 5 | trackInfo(:,5) == 3 | trackInfo(:,5) == 6));
        
        numTracksWithGap = length(unique(idxTracksWithGap));
        
        ratioGapTracks2NonGapTracks = numTracksWithGap/trackInfo(end,1);
        
        meanNNdistWithinFrame = temp.projData.medNNdistWithinFramePix.*temp.projData.pixSizeNm/1000;
        
        
        fgapFreqTime = projData.stats.fgap_freq_time;
        fgapFreqLength = projData.stats.fgap_freq_length;
        
        bgapFreqTime = projData.stats.bgap_freq_time;
        bgapFreqLength = projData.stats.bgap_freq_length;
        
        
        
        
        % find all tracks that do not contain a gap 
        % the 1/avg time or length of the growth subtrack before
        % termination of track 
        
      
           
        termFreqTime = projData.stats.term_freq_time;
        termFreqLength = projData.stats.term_freq_length;
        
        ratioFgapFreqTime2termFreqTime = fgapFreqTime/termFreqTime;
        ratioFgapFreqLength2termFreqLength =  ;
        
        ratioBgapFreqTime2termFreqTime = ;
        ratioBgapFreqLength2termFreqLength = ;
        
        ratioFgapFreqTime2BgapFreqTime = ;
        
        
        
        % Do the same for the normal per track parameters
        % assumes each project is an individual cell
        
        %MEAN VALUES
        gsAvgPerCell(iProj,1) = temp.projData.stats.growth_speed_mean_std(1,1);
        glAvgPerCell(iProj,1) = temp.projData.stats.growth_lifetime_mean_std(1,1);
        gdAvgPerCell(iProj,1) = temp.projData.stats.growth_length_mean_std(1,1);
        fsAvgPerCell(iProj,1) = temp.projData.stats.fgap_speed_mean_std(1,1);
        flAvgPerCell(iProj,1) = temp.projData.stats.fgap_lifetime_mean_std(1,1);
        fdAvgPerCell(iProj,1) = temp.projData.stats.fgap_length_mean_std(1,1);
        bsAvgPerCell(iProj,1) = temp.projData.stats.bgap_speed_mean_std(1,1);
        blAvgPerCell(iProj,1) = temp.projData.stats.bgap_lifetime_mean_std(1,1);
        bdAvgPerCell(iProj,1) = temp.projData.stats.bgap_length_mean_std(1,1);
        
      if getMedValues == 1 
        %MEDIAN VALUES
        gsMedPerCell(iProj,1) = temp.projData.stats.growth_speed_median;
        glMedPerCell(iProj,1) = temp.projData.stats.growth_lifetime_median;
        gdMedPerCell(iProj,1) = temp.projData.stats.growth_length_median;
        
        fsMedPerCell(iProj,1) = temp.projData.stats.fgap_speed_median;
        flMedPerCell(iProj,1) = temp.projData.stats.fgap_lifetime_median;
        fdMedPerCell(iProj,1) = temp.projData.stats.fgap_length_median;
        
        bsMedPerCell(iProj,1) = temp.projData.stats.bgap_speed_median;
        blMedPerCell(iProj,1) = temp.projData.stats.bgap_lifetime_median;
        bdMedPerCell(iProj,1) = temp.projData.stats.bgap_length_median;
      else 
      end % if getMedValues
      
    end
     
        % For Each Group:
        
        dataStruct(iGroup).gFor = gFor;
        dataStruct(iGroup).gBack = gBack;
        dataStruct(iGroup).gTerm = gTerm;
        dataStruct(iGroup).ratioPause2Growth = ratioPause2Growth;
        dataStruct(iGroup).ratioShrink2Growth = ratioShrink2Growth;
        dataStruct(iGroup).dynamicity = dynamicity;
        dataStruct(iGroup).ratioTotalTrack2GrowthSubtrack = ratioTotalTrack2GrowthSubtrack;
       
        % Initiate the Test ID structure so the variables can be read into
        % the discrimination matrix
        testStruct.gFor = [testID1 testID2];
        testStruct.gBack = [testID1 testID2]; 
        testStruct.gTerm = [testID1 testID2];
        testStruct.ratioPause2Growth = [testID1 testID2];
        testStruct.ratioShrink2Growth = [testID1 testID2];
        testStruct.dynamicity = [testID1 testID2];
        testStruct.ratioTotalTrack2GrowthSubtrack = [testID1 testID2];
        
        % Again do same for per track parameters (make sure to take out
        % NaN: some values particularlly bgaps will be NaN if there is no 
        % bgaps recorded in that cell)
        
        dataStruct(iGroup).gs_AvgPerCell = gsAvgPerCell(~isnan(gsAvgPerCell));
            
        dataStruct(iGroup).gl_AvgPerCell = glAvgPerCell(~isnan(glAvgPerCell));
        dataStruct(iGroup).gd_AvgPerCell = gdAvgPerCell(~isnan(gdAvgPerCell));

        dataStruct(iGroup).fs_AvgPerCell = fsAvgPerCell(~isnan(fsAvgPerCell));
        dataStruct(iGroup).fl_AvgPerCell = flAvgPerCell(~isnan(flAvgPerCell));
        dataStruct(iGroup).fd_AvgPerCell = fdAvgPerCell(~isnan(fdAvgPerCell));
        
        dataStruct(iGroup).bs_AvgPerCell = bsAvgPerCell(~isnan(bsAvgPerCell));
        dataStruct(iGroup).bl_AvgPerCell = blAvgPerCell(~isnan(blAvgPerCell));
        dataStruct(iGroup).bd_AvgPerCell = bdAvgPerCell(~isnan(bdAvgPerCell));
        
        
        
        % Per track parameters that can either be averaged per cell or
        % pooled: Agan Set up Test structures
        
        % Growth 
        testStruct.gs_AvgPerCell = [testID1 testID2];
        testStruct.gl_AvgPerCell = [testID1 testID2]; 
        testStruct.gd_AvgPerCell = [testID1 testID2];
        
        % Fgap 
        testStruct.fs_AvgPerCell = [testID1 testID2];
        testStruct.fl_AvgPerCell = [testID1 testID2];
        testStruct.fd_AvgPerCell = [testID1 testID2];
        
        % Bgap
        testStruct.bs_AvgPerCell = [testID1 testID2];
        testStruct.bl_AvgPerCell = [testID1 testID2];
           
        testStruct.bd_AvgPerCell = [testID1 testID2];
        
        if getMedValues == 1
        % Same For Median Value Parameters
        
        dataStruct(iGroup).gs_MedPerCell = gsMedPerCell(~isnan(gsMedPerCell));
        dataStruct(iGroup).gl_MedPerCell = glMedPerCell(~isnan(glMedPerCell));
        dataStruct(iGroup).gd_MedPerCell = gdMedPerCell(~isnan(gdMedPerCell));
        

        dataStruct(iGroup).fs_MedPerCell = fsMedPerCell(~isnan(fsMedPerCell));
        dataStruct(iGroup).fl_MedPerCell = flMedPerCell(~isnan(flMedPerCell));
        dataStruct(iGroup).fd_MedPerCell = fdMedPerCell(~isnan(fdMedPerCell));
        
        dataStruct(iGroup).bs_MedPerCell = bsMedPerCell(~isnan(bsMedPerCell));
        dataStruct(iGroup).bl_MedPerCell = blMedPerCell(~isnan(blMedPerCell));
        dataStruct(iGroup).bd_MedPerCell = bdMedPerCell(~isnan(bdMedPerCell));
        
        
         
        % Growth 
        testStruct.gs_MedPerCell = [testID1 testID2];
        testStruct.gl_MedPerCell = [testID1 testID2]; 
        testStruct.gd_MedPerCell = [testID1 testID2];
        
        % Fgap 
        testStruct.fs_MedPerCell = [testID1 testID2];
        testStruct.fl_MedPerCell = [testID1 testID2];
        testStruct.fd_MedPerCell = [testID1 testID2];
        
        % Bgap
        testStruct.bs_MedPerCell = [testID1 testID2];
        testStruct.bl_MedPerCell = [testID1 testID2];
        testStruct.bd_MedPerCell = [testID1 testID2];
        
        else 
        end % if getMedValues
     %%  
        
        
        % AVERAGE THE PARAMETER OVER ALL CELLS IN GROUP(for per track parameters this is an
        % average of the average value obtained per cell): 
        
        
        
        
        % parameters that can only be obtained per cell (record the mean
        % for over all cells in group)
        mean_std_gFor(iGroup,1) = mean(gFor(~isnan(gFor)));   
        mean_std_gFor(iGroup,2) = std(gFor(~isnan(gFor)));  
        
        
        mean_std_gBack(iGroup,1) = mean(gBack(~isnan(gBack)));
        mean_std_gBack(iGroup,2) = std(gBack(~isnan(gBack)));
        
        mean_std_gTerm(iGroup,1) = mean(gTerm(~isnan(gTerm)));
        
        mean_std_gTerm(iGroup,2) = std(gTerm(~isnan(gTerm)));
        
        mean_std_ratioPause2Growth(iGroup,1) = mean(ratioPause2Growth(~isnan(ratioPause2Growth)));
        mean_std_ratioPause2Growth(iGroup,2) = std(ratioPause2Growth(~isnan(ratioPause2Growth)));
        
        mean_std_ratioShrink2Growth(iGroup,1) = mean(ratioShrink2Growth(~isnan(ratioShrink2Growth)));
        mean_std_ratioShrink2Growth(iGroup,2) = std(ratioShrink2Growth(~isnan(ratioShrink2Growth)));
        
        mean_std_dynamicity(iGroup,1) = mean(dynamicity(~isnan(dynamicity)));
        mean_std_dynamicity(iGroup,2) = std(dynamicity(~isnan(dynamicity)));
        
        mean_std_ratioTotalTrack2GrowthSubtrack(iGroup,1) = mean(ratioTotalTrack2GrowthSubtrack(~isnan(ratioTotalTrack2GrowthSubtrack)));
        mean_std_ratioTotalTrack2GrowthSubtrack(iGroup,2) = std(ratioTotalTrack2GrowthSubtrack(~isnan(ratioTotalTrack2GrowthSubtrack))); 
        
        % Per track parameters that can either be averaged per cell or
        % pooled
        
        % Growth 
        mean_std_gsAvgPerCell(iGroup,1) = mean(gsAvgPerCell(~isnan(gsAvgPerCell)));
        mean_std_gsAvgPerCell(iGroup,2)= std(gsAvgPerCell(~isnan(gsAvgPerCell)));
        
        mean_std_glAvgPerCell(iGroup,1) = mean(glAvgPerCell(~isnan(glAvgPerCell)));
        mean_std_glAvgPerCell(iGroup,2) = std(glAvgPerCell(~isnan(glAvgPerCell)));
        
        mean_std_gdAvgPerCell(iGroup,1)=  mean(gdAvgPerCell(~isnan(gdAvgPerCell)));
        mean_std_gdAvgPerCell(iGroup,2) = std(gdAvgPerCell(~isnan(gdAvgPerCell)));
        
        % Fgap 
        mean_std_fsAvgPerCell(iGroup,1) = mean(fsAvgPerCell(~isnan(fsAvgPerCell)));
        mean_std_fsAvgPerCell(iGroup,2) = std(fsAvgPerCell(~isnan(fsAvgPerCell)));
        
        mean_std_flAvgPerCell(iGroup,1)= mean(flAvgPerCell(~isnan(flAvgPerCell)));
        mean_std_flAvgPerCell(iGroup,2) = std(flAvgPerCell(~isnan(flAvgPerCell)));
        
        mean_std_fdAvgPerCell(iGroup,1) = mean(fdAvgPerCell(~isnan(fdAvgPerCell)));
        mean_std_fdAvgPerCell(iGroup,2) = std(fdAvgPerCell(~isnan(fdAvgPerCell)));  
        
        % Bgap
        
        mean_std_bsAvgPerCell(iGroup,1) = mean(bsAvgPerCell(~isnan(bsAvgPerCell)));
        mean_std_bsAvgPerCell(iGroup,2) = std(bsAvgPerCell(~isnan(bsAvgPerCell)));
        
        mean_std_blAvgPerCell(iGroup,1) = mean(blAvgPerCell(~isnan(blAvgPerCell)));
        mean_std_blAvgPerCell(iGroup,2) = std(blAvgPerCell(~isnan(blAvgPerCell)));
        
        mean_std_bdAvgPerCell(iGroup,1)= mean(bdAvgPerCell(~isnan(bdAvgPerCell)));
        mean_std_bdAvgPerCell(iGroup,2) = std(bdAvgPerCell(~isnan(bdAvgPerCell)));  
        
        
        
        %% Median Values 
        
        if getMedValues == 1
        % Growth 
        mean_std_gsMedPerCell(iGroup,1) = mean(gsMedPerCell(~isnan(gsMedPerCell)));
        mean_std_gsMedPerCell(iGroup,2)= std(gsMedPerCell(~isnan(gsMedPerCell)));
        
        mean_std_glMedPerCell(iGroup,1) = mean(glMedPerCell(~isnan(glMedPerCell)));
        mean_std_glMedPerCell(iGroup,2) = std(glMedPerCell(~isnan(glMedPerCell)));
        
        mean_std_gdMedPerCell(iGroup,1)=  mean(gdMedPerCell(~isnan(gdMedPerCell)));
        mean_std_gdMedPerCell(iGroup,2) = std(gdMedPerCell(~isnan(gdMedPerCell)));
        
        % Fgap 
        mean_std_fsMedPerCell(iGroup,1) = mean(fsMedPerCell(~isnan(fsMedPerCell)));
        mean_std_fsMedPerCell(iGroup,2) = std(fsMedPerCell(~isnan(fsMedPerCell)));
        
        mean_std_flMedPerCell(iGroup,1)= mean(flMedPerCell(~isnan(flMedPerCell)));
        mean_std_flMedPerCell(iGroup,2) = std(flMedPerCell(~isnan(flMedPerCell)));
        
        mean_std_fdMedPerCell(iGroup,1) = mean(fdMedPerCell(~isnan(fdMedPerCell)));
        mean_std_fdMedPerCell(iGroup,2) = std(fdMedPerCell(~isnan(fdMedPerCell)));  
        
        % Bgap
        
        mean_std_bsMedPerCell(iGroup,1) = mean(bsMedPerCell(~isnan(bsMedPerCell)));
        mean_std_bsMedPerCell(iGroup,2) = std(bsMedPerCell(~isnan(bsMedPerCell)));
        
        mean_std_blMedPerCell(iGroup,1) = mean(blMedPerCell(~isnan(blMedPerCell)));
        mean_std_blMedPerCell(iGroup,2) = std(blMedPerCell(~isnan(blMedPerCell)));
        
        mean_std_bdMedPerCell(iGroup,1)= mean(bdMedPerCell(~isnan(bdMedPerCell)));
        mean_std_bdMedPerCell(iGroup,2) = std(bdMedPerCell(~isnan(bdMedPerCell)));  
        
        else 
        end % if getMedValues
        
        
        
end % end for iGroup

%mean_std_param is a list of the param per cell, averaged over 
% all cells (ie projects) in group
% measures the mean and the variation for the given param for all cells included 
% in the group
% Record These Values in Output Structure for Graphing

avgOverAllCellsInGroup.gFor  = mean_std_gFor;
avgOverAllCellsInGroup.gBack = mean_std_gBack;
avgOverAllCellsInGroup.gTerm = mean_std_gTerm;


avgOverAllCellsInGroup.ratioPause2Growth = mean_std_ratioPause2Growth;
avgOverAllCellsInGroup.ratioShrink2Growth = mean_std_ratioShrink2Growth;
avgOverAllCellsInGroup.dynamicity = mean_std_dynamicity;
avgOverAllCellsInGroup.ratioTotalTrack2GrowthSubtrack = mean_std_ratioTotalTrack2GrowthSubtrack;

avgOverAllCellsInGroupMean.gs = mean_std_gsAvgPerCell;
avgOverAllCellsInGroupMean.gl = mean_std_glAvgPerCell;
avgOverAllCellsInGroupMean.gd = mean_std_gdAvgPerCell;


avgOverAllCellsInGroupMean.fs = mean_std_fsAvgPerCell;
avgOverAllCellsInGroupMean.fl = mean_std_flAvgPerCell;
avgOverAllCellsInGroupMean.fd = mean_std_fdAvgPerCell;



avgOverAllCellsInGroupMean.bs = mean_std_bsAvgPerCell;
avgOverAllCellsInGroupMean.bl = mean_std_blAvgPerCell;
avgOverAllCellsInGroupMean.bd = mean_std_bdAvgPerCell;


%% Median Values

if getMedValues == 1
avgOverAllCellsInGroupMed.gs = mean_std_gsMedPerCell;
avgOverAllCellsInGroupMed.gl = mean_std_glMedPerCell;
avgOverAllCellsInGroupMed.gd = mean_std_gdMedPerCell;


avgOverAllCellsInGroupMed.fs = mean_std_fsMedPerCell;
avgOverAllCellsInGroupMed.fl = mean_std_flMedPerCell;
avgOverAllCellsInGroupMed.fd = mean_std_fdMedPerCell;



avgOverAllCellsInGroupMed.bs = mean_std_bsMedPerCell;
avgOverAllCellsInGroupMed.bl = mean_std_blMedPerCell;
avgOverAllCellsInGroupMed.bd = mean_std_bdMedPerCell;

  save([saveDir filesep 'perTrackParametersMed4Cell'],'avgOverAllCellsInGroupMed');
 
else 
end % if get MedValues

  save([saveDir filesep 'dataStruct'],'dataStruct');       
  save([saveDir filesep 'perCellParameters'],'avgOverAllCellsInGroup');
  save([saveDir filesep 'perTrackParametersMean4Cell'],'avgOverAllCellsInGroupMean');

  
  
%%
if runStats == 1
%Input the distribution of values for each project (ie cell) in group into 
% the discrimination matrix
[discrimMats]  = discriminationMatrix(dataStruct,testStruct);


pValues.gForCell = num2cell(discrimMats.gFor);
pValues.gBackCell = num2cell(discrimMats.gBack);
pValues.gTermCell = num2cell(discrimMats.gTerm);
pValues.ratioPause2GrowthCell = num2cell(discrimMats.ratioPause2Growth);
pValues.ratioShrink2GrowthCell = num2cell(discrimMats.ratioShrink2Growth);
pValues.dynamicityCell = num2cell(discrimMats.dynamicity);
pValues.ratioTotalTrack2GrowthSubtrackCell = num2cell(discrimMats.ratioTotalTrack2GrowthSubtrack);

%Per track parameters that can either be averaged per cell or
% pooled
% Growth
pValues.gs_AvgPerCell = num2cell(discrimMats.gs_AvgPerCell);
pValues.gl_AvgPerCell = num2cell(discrimMats.gl_AvgPerCell);
pValues.gd_AvgPerCell = num2cell(discrimMats.gd_AvgPerCell);

pValues.fs_AvgPerCell = num2cell(discrimMats.fs_AvgPerCell);
pValues.fl_AvgPerCell = num2cell(discrimMats.fl_AvgPerCell);
pValues.fd_AvgPerCell = num2cell(discrimMats.fd_AvgPerCell);

pValues.bs_AvgPerCell = num2cell(discrimMats.bs_AvgPerCell);
pValues.bl_AvgPerCell = num2cell(discrimMats.bl_AvgPerCell);
pValues.bd_AvgPerCell = num2cell(discrimMats.bd_AvgPerCell);





% ADD TITLES MEAN VALUES
pValues.gForCell = [groupNames pValues.gForCell];
pValues.gBackCell = [groupNames pValues.gBackCell];
pValues.gTermCell = [groupNames pValues.gTermCell];
pValues.ratioPause2GrowthCell = [groupNames pValues.ratioPause2GrowthCell];
pValues.ratioShrink2GrowthCell = [groupNames pValues.ratioShrink2GrowthCell];
pValues.dynamicityCell = [groupNames pValues.dynamicityCell];
pValues.ratioTotalTrack2GrowthSubtrackCell = [groupNames pValues.ratioTotalTrack2GrowthSubtrackCell];

pValues.gs_AvgPerCell = [groupNames pValues.gs_AvgPerCell];
pValues.gl_AvgPerCell = [groupNames pValues.gl_AvgPerCell];
pValues.gd_AvgPerCell = [groupNames pValues.gd_AvgPerCell];

pValues.fs_AvgPerCell = [groupNames pValues.fs_AvgPerCell];
pValues.fl_AvgPerCell = [groupNames pValues.fl_AvgPerCell];
pValues.fd_AvgPerCell = [groupNames pValues.fd_AvgPerCell];

pValues.bs_AvgPerCell = [groupNames pValues.bs_AvgPerCell];
pValues.bl_AvgPerCell = [groupNames pValues.bl_AvgPerCell];
pValues.bd_AvgPerCell = [groupNames pValues.bd_AvgPerCell];

if getMedValues == 1
% Median Values Per Cell
pValues.gs_MedPerCell = num2cell(discrimMats.gs_MedPerCell);
pValues.gl_MedPerCell = num2cell(discrimMats.gl_MedPerCell);
pValues.gd_MedPerCell = num2cell(discrimMats.gd_MedPerCell);

pValues.fs_MedPerCell = num2cell(discrimMats.fs_MedPerCell);
pValues.fl_MedPerCell = num2cell(discrimMats.fl_MedPerCell);
pValues.fd_MedPerCell = num2cell(discrimMats.fd_MedPerCell);

pValues.bs_MedPerCell = num2cell(discrimMats.bs_MedPerCell);
pValues.bl_MedPerCell = num2cell(discrimMats.bl_MedPerCell);
pValues.bd_MedPerCell = num2cell(discrimMats.bd_MedPerCell);




%ADD TITLES MEDIAN VALUES
pValues.gs_MedPerCell = [groupNames pValues.gs_MedPerCell];
pValues.gl_MedPerCell = [groupNames pValues.gl_MedPerCell];
pValues.gd_MedPerCell = [groupNames pValues.gd_MedPerCell];

pValues.fs_MedPerCell = [groupNames pValues.fs_MedPerCell];
pValues.fl_MedPerCell = [groupNames pValues.fl_MedPerCell];
pValues.fd_MedPerCell = [groupNames pValues.fd_MedPerCell];

pValues.bs_MedPerCell = [groupNames pValues.bs_MedPerCell];
pValues.bl_MedPerCell = [groupNames pValues.bl_MedPerCell];
pValues.bd_MedPerCell = [groupNames pValues.bd_MedPerCell];

else 
end % if getMedValues




save([saveDir filesep 'discrimMat_PerCell'],'pValues');
%%
    if getHits == 1
        
        
        %Find Hits for Second Test (first column of discrimMat)
        hitsIdx_gFor = find(discrimMats.gFor(:,1)< stringency);
        hitsIdx_gBack = find(discrimMats.gBack(:,1)<stringency);
        hitsIdx_gTerm = find(discrimMats.gTerm(:,1)<stringency);
        hitsIdx_ratioPause2Growth = find(discrimMats.ratioPause2Growth(:,1)<stringency);
        hitsIdx_ratioShrink2Growth = find(discrimMats.ratioShrink2Growth(:,1)<stringency);
        hitsIdx_dynamicity = find(discrimMats.dynamicity(:,1)<stringency);
        hitsIdx_ratioTotalTrack2GrowthSubtrack = find(discrimMats.ratioTotalTrack2GrowthSubtrack(:,1)<stringency);
        
        %Find Hits for Second Test (first column of discrimMat) 
        % PerTrackParameters
        hitsIdx_gs_AvgPerCell = find(discrimMats.gs_AvgPerCell(:,1) < stringency);
        hitsIdx_gl_AvgPerCell = find(discrimMats.gl_AvgPerCell(:,1) < stringency);
        hitsIdx_gd_AvgPerCell = find(discrimMats.gd_AvgPerCell(:,1) < stringency);
        
        hitsIdx_fs_AvgPerCell = find(discrimMats.fs_AvgPerCell(:,1) < stringency);
        hitsIdx_fl_AvgPerCell = find(discrimMats.fl_AvgPerCell(:,1) < stringency);
        hitsIdx_fd_AvgPerCell = find(discrimMats.fd_AvgPerCell(:,1) < stringency);
        
         
        hitsIdx_bs_AvgPerCell = find(discrimMats.bs_AvgPerCell(:,1) < stringency);
        hitsIdx_bl_AvgPerCell = find(discrimMats.bl_AvgPerCell(:,1) < stringency);
        hitsIdx_bd_AvgPerCell = find(discrimMats.bd_AvgPerCell(:,1) < stringency);
        
      
        
        
        
            %% Median Values: Find Hits
            if getMedValues == 1
        hitsIdx_gs_MedPerCell = find(discrimMats.gs_MedPerCell(:,1) < stringency);
        hitsIdx_gl_MedPerCell = find(discrimMats.gl_MedPerCell(:,1) < stringency);
        hitsIdx_gd_MedPerCell = find(discrimMats.gd_MedPerCell(:,1) < stringency);
        
        hitsIdx_fs_MedPerCell = find(discrimMats.fs_MedPerCell(:,1) < stringency);
        hitsIdx_fl_MedPerCell = find(discrimMats.fl_MedPerCell(:,1) < stringency);
        hitsIdx_fd_MedPerCell = find(discrimMats.fd_MedPerCell(:,1) < stringency);
        
         
        hitsIdx_bs_MedPerCell = find(discrimMats.bs_MedPerCell(:,1) < stringency);
        hitsIdx_bl_MedPerCell = find(discrimMats.bl_MedPerCell(:,1) < stringency);
        hitsIdx_bd_MedPerCell = find(discrimMats.bd_MedPerCell(:,1) < stringency);
            else 
            end 
        
      %% Generate Hit List 
       
           %Initiate gFor
           hitsList = cell(length(hitsIdx_gFor),4);
           %Write
           for iHit = 1:length(hitsIdx_gFor)
               hitsList{iHit,1} = groupNames(hitsIdx_gFor(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroup.gFor(hitsIdx_gFor(iHit),1); 
               hitsList{iHit,3} = avgOverAllCellsInGroup.gFor(hitsIdx_gFor(iHit),1)/avgOverAllCellsInGroup.gFor(1,1);
               hitsList{iHit,4} = discrimMats.gFor(hitsIdx_gFor(iHit),1);
           end
           %Save in Larger Structure
             hits.gFor= hitsList;
             
           %Initiate  gBack
           hitsList = cell(length(hitsIdx_gBack),4);
           %Write
           for iHit = 1:length(hitsIdx_gBack)
               hitsList{iHit,1} = groupNames(hitsIdx_gBack(iHit),1);
               hitsList{iHit,2} = avgOverAllCellsInGroup.gBack(hitsIdx_gBack(iHit),1);  
               hitsList{iHit,3} = avgOverAllCellsInGroup.gBack(hitsIdx_gBack(iHit),1)/avgOverAllCellsInGroup.gBack(1,1);
               hitsList{iHit,4} = discrimMats.gBack(hitsIdx_gBack(iHit),1);
           end
           %Save in Larger Structure
           hits.gBack = hitsList;
                
           %Initiate gTerm
           hitsList = cell(length(hitsIdx_gTerm),4);
           % Write
           for iHit = 1:length(hitsIdx_gTerm)
               hitsList{iHit,1} = groupNames(hitsIdx_gTerm(iHit),1);
               hitsList{iHit,2} = avgOverAllCellsInGroup.gTerm(hitsIdx_gTerm(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroup.gTerm(hitsIdx_gTerm(iHit),1)/avgOverAllCellsInGroup.gTerm(1,1);
               hitsList{iHit,4} = discrimMats.gTerm(hitsIdx_gTerm(iHit),1);
           end
           %Save in Larger Structure
           hits.gTerm = hitsList;
           
           
           %Initiate Pause2Growth
           hitsList = cell(length(hitsIdx_ratioPause2Growth),4);
           %Write
           for iHit = 1:length(hitsIdx_ratioPause2Growth)
              hitsList{iHit,1} = groupNames(hitsIdx_ratioPause2Growth(iHit),1);
              hitsList{iHit,2} = avgOverAllCellsInGroup.ratioPause2Growth(hitsIdx_ratioPause2Growth(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroup.ratioPause2Growth(hitsIdx_ratioPause2Growth(iHit),1)/avgOverAllCellsInGroup.ratioPause2Growth(1,1);
              hitsList{iHit,4} = discrimMats.ratioPause2Growth(hitsIdx_ratioPause2Growth(iHit),1);           
           end 
           
           %Save in Larger Structure
           hits.ratioPause2Growth = hitsList;
           
           %Initiate ratioShrink2Growth 
           hitsList = cell(length(hitsIdx_ratioShrink2Growth),4);
           %Write
           for iHit = 1:length(hitsIdx_ratioShrink2Growth) 
               hitsList{iHit,1} = groupNames(hitsIdx_ratioShrink2Growth(iHit),1);
               hitsList{iHit,2} = avgOverAllCellsInGroup.ratioShrink2Growth(hitsIdx_ratioShrink2Growth(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroup.ratioShrink2Growth(hitsIdx_ratioShrink2Growth(iHit),1)/avgOverAllCellsInGroup.ratioShrink2Growth(1,1);
               hitsList{iHit,4} = discrimMats.ratioShrink2Growth(hitsIdx_ratioShrink2Growth(iHit),1);
           end 
           %Save in Larger Structure
              hits.ratioShrink2Growth = hitsList;
           
           %Initiate dynamicity 
           hitsList = cell(length(hitsIdx_dynamicity),4);
           %Write
           for iHit = 1:length(hitsIdx_dynamicity) 
               hitsList{iHit,1} = groupNames(hitsIdx_dynamicity(iHit),1);
               hitsList{iHit,2} = avgOverAllCellsInGroup.dynamicity(hitsIdx_dynamicity(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroup.dynamicity(hitsIdx_dynamicity(iHit),1)/avgOverAllCellsInGroup.dynamicity(1,1);
               hitsList{iHit,4} = discrimMats.dynamicity(hitsIdx_dynamicity(iHit),1);
           end 
           %Save in Larger Structure
              hits.dynamicity = hitsList;
              
           %Initiate ratioTotalTrack2GrowthSubtrack
           hitsList = cell(length(hitsIdx_ratioTotalTrack2GrowthSubtrack),4);
           %Write
           for iHit = 1:length(hitsIdx_ratioTotalTrack2GrowthSubtrack) 
               hitsList{iHit,1} = groupNames(hitsIdx_ratioTotalTrack2GrowthSubtrack(iHit),1);
               hitsList{iHit,2} = avgOverAllCellsInGroup.ratioTotalTrack2GrowthSubtrack(hitsIdx_ratioTotalTrack2GrowthSubtrack(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroup.ratioTotalTrack2GrowthSubtrack(hitsIdx_ratioTotalTrack2GrowthSubtrack(iHit),1)/avgOverAllCellsInGroup.ratioTotalTrack2GrowthSubtrack(1,1);
               hitsList{iHit,4} = discrimMats.ratioTotalTrack2GrowthSubtrack(hitsIdx_ratioTotalTrack2GrowthSubtrack(iHit),1);
           end 
           %Save in Larger Structure
              hits.ratioCompound2Growth = hitsList; 
              
              
              
          %% Generate Hits List: Per Track Parameters- Avg Per Cell
           % Initiate gs_AvgPerCell
           hitsList = cell(length(hitsIdx_gs_AvgPerCell),4);
           %Write
           for iHit = 1:length(hitsIdx_gs_AvgPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_gs_AvgPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMean.gs(hitsIdx_gs_AvgPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMean.gs(hitsIdx_gs_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.gs(1,1);
              hitsList{iHit,4} = discrimMats.gs_AvgPerCell(hitsIdx_gs_AvgPerCell(iHit),1);
           end
           %Save to larger structure
           hits.gs_AvgPerCell = hitsList;
           
           
           % Initiate gl_AvgPerCell
           hitsList = cell(length(hitsIdx_gl_AvgPerCell),4);
           %Write
           for iHit = 1:length(hitsIdx_gl_AvgPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_gl_AvgPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMean.gl(hitsIdx_gl_AvgPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMean.gl(hitsIdx_gl_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.gl(1,1);
               hitsList{iHit,4} = discrimMats.gl_AvgPerCell(hitsIdx_gl_AvgPerCell(iHit),1);
               
           end
           hits.gl_AvgPerCell = hitsList;
           
           % Initiate gd_AvgPerCell
           hitsList = cell(length(hitsIdx_gd_AvgPerCell),4);
           
           %Write
           for iHit = 1:length(hitsIdx_gd_AvgPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_gd_AvgPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMean.gd(hitsIdx_gd_AvgPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMean.gd(hitsIdx_gd_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.gd(1,1);
               hitsList{iHit,4} = discrimMats.gd_AvgPerCell(hitsIdx_gd_AvgPerCell(iHit),1);
           end
           
            hits.gd_AvgPerCell = hitsList;
           
           
           
           
           %Fgap parameters
           hitsList = cell(length(hitsIdx_fs_AvgPerCell),4);
           
            for iHit = 1:length(hitsIdx_fs_AvgPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_fs_AvgPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMean.fs(hitsIdx_fs_AvgPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMean.fs(hitsIdx_fs_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.fs(1,1);
              hitsList{iHit,4} = discrimMats.fs_AvgPerCell(hitsIdx_fs_AvgPerCell(iHit),1);
            end
           
            hits.fs_AvgPerCell = hitsList;
            
            % Initiate fl_AvgPerCell
            hitsList = cell(length(hitsIdx_fl_AvgPerCell),4);
           for iHit = 1:length(hitsIdx_fl_AvgPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_fl_AvgPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMean.fl(hitsIdx_fl_AvgPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMean.fl(hitsIdx_fl_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.fl(1,1);
               hitsList{iHit,4} = discrimMats.fl_AvgPerCell(hitsIdx_fl_AvgPerCell(iHit),1);
           end
           hits.fl_AvgPerCell = hitsList;
           
           % Initiate fd_AvgPerCell
            hitsList = cell(length(hitsIdx_fd_AvgPerCell),4);
            
           for iHit = 1:length(hitsIdx_fd_AvgPerCell)
             hitsList{iHit,1} = groupNames(hitsIdx_fd_AvgPerCell(iHit),1); 
             hitsList{iHit,2} = avgOverAllCellsInGroupMean.fd(hitsIdx_fd_AvgPerCell(iHit),1);
             hitsList{iHit,3} = avgOverAllCellsInGroupMean.fd(hitsIdx_fd_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.fd(1,1);
             hitsList{iHit,4} = discrimMats.fd_AvgPerCell(hitsIdx_fd_AvgPerCell(iHit),1);
           end
              hits.fd_AvgPerCell = hitsList;
            
           %Bgap parameters
           
           hitsList = cell(length(hitsIdx_bs_AvgPerCell),4);
           
            for iHit = 1:length(hitsIdx_bs_AvgPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_bs_AvgPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMean.bs(hitsIdx_bs_AvgPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMean.bs(hitsIdx_bs_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.bs(1,1);
              hitsList{iHit,4} = discrimMats.bs_AvgPerCell(hitsIdx_bs_AvgPerCell(iHit),1);
           end
           hits.bs_AvgPerCell = hitsList;
           
           hitsList = cell(length(hitsIdx_bl_AvgPerCell),4);
           for iHit = 1:length(hitsIdx_bl_AvgPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_bl_AvgPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMean.bl(hitsIdx_bl_AvgPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMean.bl(hitsIdx_bl_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.bl(1,1);
               hitsList{iHit,4} = discrimMats.bl_AvgPerCell(hitsIdx_bl_AvgPerCell(iHit),1);
           end
           hits.bl_AvgPerCell = hitsList;
           
           hitsList = cell(length(hitsIdx_bd_AvgPerCell),4);
           for iHit = 1:length(hitsIdx_bd_AvgPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_bd_AvgPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMean.bd(hitsIdx_bd_AvgPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMean.bd(hitsIdx_bd_AvgPerCell(iHit),1)/avgOverAllCellsInGroupMean.bd(1,1);
               hitsList{iHit,4} = discrimMats.bd_AvgPerCell(hitsIdx_bd_AvgPerCell(iHit),1);
           end
           hits.bd_AvgPerCell = hitsList;
           
           
           
           %% Generate Hit List: Median Values 
           if getMedValues == 1
           %Growth Parameters
           % Initiate gs_MedPerCell
           hitsList = cell(length(hitsIdx_gs_MedPerCell),4);
           %Write
           for iHit = 1:length(hitsIdx_gs_MedPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_gs_MedPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMed.gs(hitsIdx_gs_MedPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMed.gs(hitsIdx_gs_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.gs(1,1);
              hitsList{iHit,4} = discrimMats.gs_MedPerCell(hitsIdx_gs_MedPerCell(iHit),1);
           end
           %Save to larger structure
           hitsMed.gs_MedPerCell = hitsList;
           
           
           % Initiate gl_AvgPerCell
           hitsList = cell(length(hitsIdx_gl_MedPerCell),4);
           %Write
           for iHit = 1:length(hitsIdx_gl_MedPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_gl_MedPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMed.gl(hitsIdx_gl_MedPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMed.gl(hitsIdx_gl_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.gl(1,1);
               hitsList{iHit,4} = discrimMats.gl_MedPerCell(hitsIdx_gl_MedPerCell(iHit),1);
               
           end
           hitsMed.gl_MedPerCell = hitsList;
           
           % Initiate gd_AvgPerCell
           hitsList = cell(length(hitsIdx_gd_MedPerCell),4);
           
           %Write
           for iHit = 1:length(hitsIdx_gd_MedPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_gd_MedPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMed.gd(hitsIdx_gd_MedPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMed.gd(hitsIdx_gd_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.gd(1,1);
               hitsList{iHit,4} = discrimMats.gd_MedPerCell(hitsIdx_gd_MedPerCell(iHit),1);
           end
           
            hitsMed.gd_MedPerCell = hitsList;
           
           
           
           
           %Fgap parameters
           hitsList = cell(length(hitsIdx_fs_MedPerCell),4);
           
            for iHit = 1:length(hitsIdx_fs_MedPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_fs_MedPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMed.fs(hitsIdx_fs_MedPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMed.fs(hitsIdx_fs_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.fs(1,1);
              hitsList{iHit,4} = discrimMats.fs_AvgPerCell(hitsIdx_fs_MedPerCell(iHit),1);
            end
           
            hitsMed.fs_MedPerCell = hitsList;
            
            % Initiate fl_AvgPerCell
            hitsList = cell(length(hitsIdx_fl_MedPerCell),4);
           for iHit = 1:length(hitsIdx_fl_MedPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_fl_MedPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMed.fl(hitsIdx_fl_MedPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMed.fl(hitsIdx_fl_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.fl(1,1);
               hitsList{iHit,4} = discrimMats.fl_AvgPerCell(hitsIdx_fl_MedPerCell(iHit),1);
           end
           hitsMed.fl_MedPerCell = hitsList;
           
           % Initiate fd_AvgPerCell
            hitsList = cell(length(hitsIdx_fd_MedPerCell),4);
            
           for iHit = 1:length(hitsIdx_fd_MedPerCell)
             hitsList{iHit,1} = groupNames(hitsIdx_fd_MedPerCell(iHit),1); 
             hitsList{iHit,2} = avgOverAllCellsInGroupMed.fd(hitsIdx_fd_MedPerCell(iHit),1);
             hitsList{iHit,3} = avgOverAllCellsInGroupMed.fd(hitsIdx_fd_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.fd(1,1);
             hitsList{iHit,4} = discrimMats.fd_MedPerCell(hitsIdx_fd_MedPerCell(iHit),1);
           end
              hitsMed.fd_MedPerCell = hitsList;
            
           %Bgap parameters
           
           hitsList = cell(length(hitsIdx_bs_MedPerCell),4);
           
            for iHit = 1:length(hitsIdx_bs_MedPerCell)
              hitsList{iHit,1} = groupNames(hitsIdx_bs_MedPerCell(iHit),1); 
              hitsList{iHit,2} = avgOverAllCellsInGroupMed.bs(hitsIdx_bs_MedPerCell(iHit),1);
              hitsList{iHit,3} = avgOverAllCellsInGroupMed.bs(hitsIdx_bs_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.bs(1,1);
              hitsList{iHit,4} = discrimMats.bs_MedPerCell(hitsIdx_bs_MedPerCell(iHit),1);
           end
           hitsMed.bs_MedPerCell = hitsList;
           
           hitsList = cell(length(hitsIdx_bl_MedPerCell),4);
           for iHit = 1:length(hitsIdx_bl_MedPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_bl_MedPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMed.bl(hitsIdx_bl_MedPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMed.bl(hitsIdx_bl_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.bl(1,1);
               hitsList{iHit,4} = discrimMats.bl_MedPerCell(hitsIdx_bl_MedPerCell(iHit),1);
           end
           hitsMed.bl_MedPerCell = hitsList;
           
           hitsList = cell(length(hitsIdx_bd_MedPerCell),4);
           for iHit = 1:length(hitsIdx_bd_MedPerCell)
               hitsList{iHit,1} = groupNames(hitsIdx_bd_MedPerCell(iHit),1); 
               hitsList{iHit,2} = avgOverAllCellsInGroupMed.bd(hitsIdx_bd_MedPerCell(iHit),1);
               hitsList{iHit,3} = avgOverAllCellsInGroupMed.bd(hitsIdx_bd_MedPerCell(iHit),1)/avgOverAllCellsInGroupMed.bd(1,1);
               hitsList{iHit,4} = discrimMats.bd_MedPerCell(hitsIdx_bd_MedPerCell(iHit),1);
           end
           hitsMed.bd_MedPerCell = hitsList;
           save([saveDir filesep 'hitsMed'],'hitsMed');
           
           else 
           end % if getMedValues
           
           
           
           
           
           
    save([saveDir filesep 'hits'],'hits');
    
    
   
    end % end if get hits
end % end if runStats
        
   
        

end 



