function targetStats(mData)

doPlot=1;

% get output directory and excel sheet for writing
homeDir=pwd;
outputDir=uigetdir(pwd,'Select output directory.');
cd(outputDir)
filePath=[];
while exist(filePath,'file')==0
    [fileName,pathName]=uigetfile('*.xls','Select .xls file for data output.');
    filePath=[pathName fileName];
end
cd(homeDir)

% use doPlot variable for output directory path
if doPlot==1
    doPlot=outputDir;
else
    doPlot=[];
end

lastRow=0;

% get list of all target names
targetLabels=getlabels(mData.target);


% exList=[];
% for iTar=1:2
%     
%     [movGroup]=splitByMovie(mData,iTar,[]);
% 
%     discrimMatrices=[];
%     prop2sample='growthSpeeds';
%     [movGroup]=getPooledData(mData,movGroup,prop2sample);
%     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);
% 
% %     prop2sample='shrinkSpeeds';
% %     [movGroup]=getPooledData(mData,movGroup,prop2sample);
% %     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);
% 
%     [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow);
% 
%     save([outputDir filesep 'targetStatsOutput1 ' num2str(iTar)],'movGroup','discrimMatrices')
% 
% end
lastRow=134;
exList=46;

for iTar=1:2

    [movGroup]=splitByDate(mData,iTar,exList);

    discrimMatrices=[];
    prop2sample='growthSpeeds';
    [movGroup]=getPooledData(mData,movGroup,prop2sample);
    [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

%     prop2sample='shrinkSpeeds';
%     [movGroup]=getPooledData(mData,movGroup,prop2sample);
%     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

    [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow);

    save([outputDir filesep 'targetStatsOutput2 ' num2str(iTar)],'movGroup','discrimMatrices')

end

for iTar=3:length(targetLabels);
    [movGroup]=splitByOligo(mData,iTar,exList);

    discrimMatrices=[];
    prop2sample='growthSpeeds';
    [movGroup]=getPooledData(mData,movGroup,prop2sample);
    [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

%     prop2sample='shrinkSpeeds';
%     [movGroup]=getPooledData(mData,movGroup,prop2sample);
%     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);


    [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow);

    save([outputDir filesep 'targetStatsOutput3 ' num2str(iTar)],'movGroup','discrimMatrices')

end




%%%%% sub-functions %%%%%
function [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow)

letters={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'...
    'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'...
    'AA' 'AB' 'AC' 'AD' 'AE' 'AF' 'AG' 'AH' 'AI' 'AJ' 'AK' 'AL' 'AM' 'AN'...
    'AO' 'AP' 'AQ' 'AR' 'AS' 'AT' 'AU' 'AV' 'AW' 'AX' 'AY' 'AZ'...
    'BA' 'BB' 'BC' 'BD' 'BE' 'BF' 'BG' 'BH' 'BI' 'BJ' 'BK' 'BL' 'BM' 'BN'...
    'BO' 'BP' 'BQ' 'BR' 'BS' 'BT' 'BU' 'BV' 'BW' 'BX' 'BY' 'BZ'...
    'CA' 'CB' 'CC' 'CD' 'CE' 'CF' 'CG' 'CH' 'CI' 'CJ' 'CK' 'CL' 'CM' 'CN'...
    'CO' 'CP' 'CQ' 'CR' 'CS' 'CT' 'CU' 'CV' 'CW' 'CX' 'CY' 'CZ'};

names = fieldnames(discrimMatrices);
for iName=1:length(names)
    tempMat=discrimMatrices.(names{iName});
    if iscell(tempMat)
        [r c]=size(tempMat);
        range=[letters{2} num2str(lastRow+3) ':' letters{1+c} num2str(lastRow+2+r)];
        xlswrite(filePath,tempMat,range);
        lastRow=lastRow+2+r;
    end

end


function [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot)

switch prop2sample
    case {'growthSpeeds','shrinkSpeeds'} % average KS (mean subtracted) test on pop, t-test on sample means

        nReps=100; % do 100 repetitions
        maxSampleSize=400; % max sample size is 400 for KS test
        nGroups=length(movGroup); % number of groups to sample

        % nGroups-vector of population sizes
        popSize=arrayfun(@(x)length(x.(prop2sample)),movGroup);

        % sample each group with smallest distribution size
        sampSizeLargest=min(popSize);
        sampleDataLargest=arrayfun(@(x)randsample(x.(prop2sample),sampSizeLargest),movGroup,'UniformOutput', false);

        if ~isempty(doPlot)
            % make boxplot of group distributions
            figure
            boxplot(cell2mat(sampleDataLargest'))
            xlabel(' ');
            ylabel(prop2sample);
            set(gca,'XTickLabel',{movGroup.label}');
            rotateticklabel(gca,30);
            title([movGroup(1,1).common2group ' - ' prop2sample]);
            
            plotDir=[doPlot filesep 'boxPlots'];
            if ~isdir(plotDir)
                mkdir(plotDir)
            end
            saveas(gcf,[plotDir filesep movGroup(1,1).common2group '_' prop2sample '.fig']);
            close all
        end

        % sample each group for maxSampleSize values
        sampSizeKS=min([sampSizeLargest; maxSampleSize]);

        %initialize discrimination matrix
        tempDiscrimMat=zeros(nGroups,nGroups,nReps);

        % KS test with mean subtracted = 11
        testStructure.(prop2sample)=[11 11];

        sMeans=zeros(nReps,nGroups);
        for iRep=1:nReps
            % get maxSampleSize values from each group, put in structure, and find
            % mean to be used in later test
            sampleDataKS=arrayfun(@(x)randsample(x.(prop2sample),sampSizeKS),movGroup,'UniformOutput', false);
            sData=cell2struct(sampleDataKS,prop2sample,2);
            sMeans(iRep,:)=mean([sData.(prop2sample)]);

            % make discrimination matrix for property - KS test of sampled population
            compMatrices=discriminationMatrix(sData,testStructure);
            tempDiscrimMat(:,:,iRep)=compMatrices.(prop2sample);
        end

        M=diag(ones(nGroups,1)); M=swapMaskValues(M,[0 1],[1 nan]);
        lowerTri=tril(ones(nGroups));
        upperTri=triu(ones(nGroups));

        % get average discrimination matrix for property from nRep trials
        discrimMatValues=mean(tempDiscrimMat,3).*M;

        clear testStructure
        clear sData

        % make a structure containing mean values from the nReps trials
        sData=cell2struct(mat2cell(sMeans,nReps,ones(nGroups,1))',[prop2sample 'Means'],2);

        % make discrimination matrix for property - t-test of sampled means
        testStructure.([prop2sample 'Means'])=[1 1];
        tempDiscrimMat=discriminationMatrix(sData,testStructure);
        discrimMatMeans=tempDiscrimMat.([prop2sample 'Means']).*M;

        % store disc matrices in structure
        discrimMatrices.(prop2sample)=discrimMatValues;
        discrimMatrices.([prop2sample 'Means'])=discrimMatMeans;
        
        % combine the two matrices into one, with labels
        tempMerge=num2cell(discrimMatValues.*upperTri + discrimMatMeans.*lowerTri);
        idx=find(eye(size(tempMerge,1)));
        tempMerge(idx)=cellfun(@(x)num2str(x),tempMerge(idx),'UniformOutput',false);
        
        cLabel={movGroup.label};
        rLabel=[{prop2sample}; {movGroup.label}'];
        tempMerge=[rLabel [cLabel; tempMerge]];
        discrimMatrices.([prop2sample 'TKS'])=tempMerge;
        
    otherwise
        error('prop2sample not supported')
end



function [movGroup]=splitByOligo(mData,iTar,exList)
if nargin<3
    exList=[];
end

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

% get list of all oligo names corresponding to iTar
iTarOligoLabels=getlabels(nominal(mData.oligo(movNums)));

% loop thru oligo names and find corresponding movie numbers, and remove
% any that are to be excluded
nOligos=length(iTarOligoLabels);
for iOligo=1:nOligos
    movGroup(iOligo,1).common2group=[targetLabels{iTar} '_split_by_oligo']; 
    movGroup(iOligo,1).label=[targetLabels{iTar} '_' iTarOligoLabels{iOligo}];
    % movies are split by oligo number only - same oligos on different days
    % will be grouped together
    movNums=find(ismember(mData.oligo,iTarOligoLabels(iOligo)))';
    movGroup(iOligo,1).movNums=setdiff(movNums,exList)';
end


function [movGroup]=splitByMovie(mData,iTar,exList)
if nargin<3
    exList=[];
end

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

for iMov=1:length(movNums)
    movGroup(iMov,1).common2group=[targetLabels{iTar} '_split_by_movie'];
    movGroup(iMov,1).label=[targetLabels{iTar} '_' mData.oligo{movNums(iMov)}...
        '_' mData.date{movNums(iMov)} '_mov_' num2str(movNums(iMov))];
    movGroup(iMov,1).movNums=movNums(iMov);
       
end





% this function mainly good for grouping control movies by date
function [movGroup]=splitByDate(mData,iTar,exList)
if nargin<3
    exList=[];
end

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

% get list of all date names corresponding to iTar
iTarDateLabels=getlabels(nominal(mData.date(movNums)));

% loop thru date names and find corresponding movie numbers, and remove
% any that are to be excluded
nDates=length(iTarDateLabels);
for iDate=1:nDates
    movGroup(iDate,1).common2group=[targetLabels{iTar} '_split_by_date'];
    movGroup(iDate,1).label=[targetLabels{iTar} '_' iTarDateLabels{iDate}];
    % movies must be split by date and belong to target
    movNums=find(ismember(mData.date,iTarDateLabels(iDate)) & ismember(mData.target,targetLabels(iTar)))';
    movGroup(iDate,1).movNums=setdiff(movNums,exList)';
end

function [movGroup]=getPooledData(mData,movGroup,prop2sample)

nGroups=length(movGroup);
for iGroup=1:nGroups
    nMovs=length(movGroup(iGroup,1).movNums);
    pooledData=cell(nMovs,1);
    for iMov=1:nMovs
        % load projData for iMovie
        dataLoc=[mData{movGroup(iGroup,1).movNums(iMov,1),6} filesep 'meta'];
        load([dataLoc filesep 'projData']);

        switch prop2sample
            case 'growthSpeeds'
                % extract the track profile matrix
                tempProf=projData.nTrack_start_end_velMicPerMin_class_lifetime;
                % extract growth speeds
                iMovData=tempProf(tempProf(:,5)==1,4);
            case 'shrinkSpeeds'
                % extract the track profile matrix
                tempProf=projData.nTrack_start_end_velMicPerMin_class_lifetime;
                % extract shrinkage speeds
                iMovData=tempProf(tempProf(:,5)==3,4);
            otherwise
                error('prop2sample not supported')
        end
        pooledData{iMov,1}=iMovData;
    end
    movGroup(iGroup,1).(prop2sample)=cell2mat(pooledData);
end



% g1= tempDataset.trackID<=10 & tempDataset.trackID<20;
% find(g1)
% bigSet=vertcat(trackProfiles{1:3,1});
% getlabels(bigSet.movID)
% group=bigSet.trackType
% [order,number,group_median,group_iqr] = grpstats(bigSet.vel,group,{'gname','numel',@median,@iqr})
%
%

% nTracks=size(tempProf,1);
%
%         movID = nominal(repmat(movGroup.movNums(iMov),[nTracks,1]));
% %         date  = ordinal(repmat(mData{movGroup.movNums(iMov),1},[nTracks,1]));
% %         target=        (repmat(mData{movGroup.movNums(iMov),2},[nTracks,1])); % already nominal
% %         oligo = nominal(repmat(mData{movGroup.movNums(iMov),3},[nTracks,1]));
% %         movNum= nominal(repmat(mData{movGroup.movNums(iMov),4},[nTracks,1]));
% %         roiNum= nominal(repmat(mData{movGroup.movNums(iMov),5},[nTracks,1]));
%         %{date,'date'},{target,'target'},{oligo,'oligo'},{movNum,'movNum'},{roiNum,'roiNum'},
%         tempDataset=dataset({movID,'movID'},...
%             {tempProf,'trackID','sFr','eFr','vel','trackType','length'});
%
%         unique(tempDataset.trackType);
%
%         tempDataset.trackType=nominal(tempDataset.trackType,...
%             {'NA','growth','fgap','shrink','unclass'},[],[-.5:1:4.5]);
%
%         trackProfiles{iMov,1}=tempDataset;

