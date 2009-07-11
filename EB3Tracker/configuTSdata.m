function [allData]=configuTSdata(mData,exList,allData)
%close all
% trackNum=zeros(100,1);
% for i=1:100
%
%     temp=load([formatPath(projList(i,1).anDir) filesep 'meta' filesep 'projData']);
%     projData=temp.projData;
%
%     trackNum(i)=projData.numTracks;
% end
%
%
% % first, load projList and excludeListDensity, which was saved by taking
% % 99% of stds from the plusTipCheckFeatureDensity function.
% [mData]=generateMovieDatabase(projList);
% targetLabels=getlabels(mData.target);
% [allProjects,notDone]=plusTipCheckIfDone(2);
% %Selected: V:\EB3_TS_Screen\activeData\\projListThru090708.mat
% exList=[excludeListDensity; notDone; 46]; % 46 is noticeably off from the others Allstars



% % movNumsVHL=find(mData.target=='VHL');
% % movNumsAllstar=[1:4]';
% % movNumsNegctr=[17:20]';
%
% movNums=find(mData.target=='Negctr'); %[movNumsVHL; movNumsAllstar; movNumsNegctr];
%
%
% for i=1:20 %length(movNums)
%     temp=load([formatPath(mData{movNums(i),7}) filesep 'meta' filesep 'projData']);
%     projData=temp.projData;
%
%     data=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
%
%     growthSpeed{i,1}=data(data(:,5)==1,4);
%     growthLFT{i,1}=data(data(:,5)==1,6);
%     shrinkSpeed{i,1}=abs(data(data(:,5)==3,4));
%     medNNnm{i,1}=projData.medNNdistWithinFramePix.*projData.pixSizeNm./1000;
%
% end
%
% minVal=min(cellfun(@(x) min(x),growthSpeed));
% maxVal=max(cellfun(@(x) prctile(x,95),growthSpeed));
%
% nTracks=cellfun(@(x) length(x),growthSpeed);
%
% n=linspace(minVal,maxVal,40);
%
% [x1,binList]=cellfun(@(x,y) histc(x,n),growthSpeed,'uniformOutput',0);
%
% t=cell2mat(x1');
% tnorm=t./repmat(nTracks',[40 1]);
% cMap=[repmat([1 0 0],12,1); repmat([0 1 0],4,1); repmat([0 0 1],4,1)];
%
% figure; h=bar3(n,tnorm);
% set(h,'CData',repmat((1:6*size(t,1)).',1,4));



%


% targetLabels=getlabels(mData.target);
% for iTar=1:15 %length(targetLabels);
%     [movGroup]=splitByOligo(mData,iTar,exList);
%     allData(iTar,1).movGroup=movGroup;
%
% end
%
% for iTar=1:length(allData)
%     nOligo=length(allData(iTar,1).movGroup);
%     for iOligo=1:nOligo
%         movs=allData(iTar,1).movGroup(iOligo,1).movNums;
%         for i=1:length(movs)
%             temp=load([formatPath(mData{movs(i),7}) filesep 'meta' filesep 'projData']);
%             projData=temp.projData;
%
%             data=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
%
%             allData(iTar,1).movGroup(iOligo,1).growthSpeed{i,1}=data(data(:,5)==1,4);
%             allData(iTar,1).movGroup(iOligo,1).growthLFT{i,1}=data(data(:,5)==1,6);
%             allData(iTar,1).movGroup(iOligo,1).pauseSpeed{i,1}=abs(data(data(:,5)==2,4));
%             allData(iTar,1).movGroup(iOligo,1).shrinkSpeed{i,1}=abs(data(data(:,5)==3,4));
%             allData(iTar,1).movGroup(iOligo,1).medNNnm{i,1}=projData.medNNdistWithinFramePix.*projData.pixSizeNm./1000;
%
%         end
%     end
% end
% save('allData','allData')



paramList={'growthSpeed','pauseSpeed','shrinkSpeed'};
for iParam=1:length(paramList)
    count=1;
    for tar=1:14
        % count how many oligos exist per gene
        nOli=size(allData(tar,1).movGroup,1);
        % get all data values in an nOligo-long cell array
        oli=arrayfun(@(t,m) cell2mat(allData(t,1).movGroup(m,1).(paramList{iParam})),tar.*ones(1,nOli),1:nOli,'uniformoutput',0)';

        % figure out how many values to use for sampling (need same number for
        % each condition, so use the smallest sample)
        if tar~=2
            sampleSize=min([10000; cell2mat(cellfun(@(x) length(x),oli,'uniformoutput',0))]);
        else
            sampleSize=3*sampleSize; % negctrl only has one oligo
        end


        tempVal=cellfun(@(x) randsample(x,sampleSize),oli,'uniformoutput',0)';
        tempLab={allData(tar,1).movGroup(:,1).label};

        oliVal(1,count:count+nOli-1)=tempVal;
        oliLab(1,count:count+nOli-1)=tempLab;

        tarAll{tar,1}=reshape(cell2mat(tempVal),[],1);
        tarLab{tar,1}=allData(tar,1).movGroup.common2group;

        count=count+nOli;
    end

    % plot the data by oligo
    % figure; boxplot(cell2mat(oliVal),'notch','on')
    % set(gca,'XTickLabel',oliLab);
    % rotateticklabel(gca,45);

    % cut down total N to smallest number represented by any of the
    % conditions
    nVal=cellfun(@(x) length(x),tarAll);
    idx=find(nVal>min(nVal));
    for j=1:length(idx)
        tarAll{idx(j),1}(min(nVal)+1:end)=[];
    end

    % plot the data by target
    %figure
    subplot(3,1,iParam)
    boxplot(cell2mat(tarAll'),'notch','on')
    set(gca,'XTickLabel',tarLab');
    
    switch iParam
        case 1
            ylabel('Growth Speed')
            ylim([0 40])
%         case 2
%             ylabel('Growth Phase Lifetime (frames)')
%             ylim([0 25])
        case 2
            ylabel('Pause Speed')
            ylim([0 25])
        case 3
            ylabel('Shrink Speed')
            ylim([0 20])            
    end

end




tarAll;




function [movGroup]=splitByOligo(mData,iTar,exList)
if nargin<3
    exList=[];
end
movGroup=[];

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
    movGroup(iOligo,1).common2group=[targetLabels{iTar}];
    movGroup(iOligo,1).label=[targetLabels{iTar} '_' iTarOligoLabels{iOligo}];
    % movies are split by oligo number only - same oligos on different days
    % will be grouped together
    movNums=find(ismember(mData.oligo,iTarOligoLabels(iOligo)))';
    movGroup(iOligo,1).movNums=setdiff(movNums,exList)';
end




