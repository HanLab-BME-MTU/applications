function plusTipStackPercentBars(groupList)


dirName=uigetdir(pwd,'Choose parent directory containing all group quad plots');

if ispc
    fileExt='.emf';
else
    fileExt='.jpg';
end

% get sorted group names in order of appearance
[groupNames,m,n]=unique(groupList(:,1));
idx=sort(m);
groupNames=groupList(idx,1);

%groupNames=cellfun(@(x) ['grp_' x],groupNames,'uniformoutput',0);

% put populations into a cell array
for i=1:length(groupNames)
fileName=[dirName filesep groupNames{i} '_groupQuadInfo.mat'];
    if exist(fileName,'file')~=0
       load([dirName filesep groupNames{i} '_groupQuadInfo'])
    else
        error(['No such file as' fileName]);        
    end
   allPop{i,1}=popRYGB;
end
% sum the total tracks in each group
sumPopPerGroup=cell2mat(cellfun(@(x) sum(x),allPop,'uniformOutput',0));
% get the percent bar
[nPrctRYGB]=plusTipQuadColorbar(sumPopPerGroup);

save([dirName filesep 'stackedPrctBarNumbers'],'groupNames','allPop','sumPopPerGroup','nPrctRYGB')
saveas(gcf,[dirName filesep 'prctBarStackedGroups' '.fig'])
saveas(gcf,[dirName filesep 'prctBarStackedGroups' fileExt])
