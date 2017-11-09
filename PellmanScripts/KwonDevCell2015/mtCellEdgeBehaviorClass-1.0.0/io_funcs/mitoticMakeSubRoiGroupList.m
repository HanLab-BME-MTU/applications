function [ groupList] = mitoticMakeSubRoiGroupList( groupList ,varargin)
% mitoticMakeSubRoiGroupList.m 
% 
% 
%



%% CheckInput
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true ; 
ip.addParameter('subRoiFilename',[], @(x) ischar(x) || isempty(x));

ip.addParameter('bipolar',false,@(x) islogical(x)); % in um

ip.addParameter('sortGroupList',true, @(x) islogical(x));

ip.addParameter('outputDirectory',[]); 

ip.addParameter('subRoiGroupListFilename',[]); 

ip.addParameter('replaceStr',[]); 

ip.parse(varargin{:});
%% Set Up 
% groupList filename 

if isempty(ip.Results.subRoiFilename)
    if ip.Results.bipolar
        subRoiFilename = [ 'mitoticSubRois_1_um_bipolar'];
    else
        subRoiFilename = ['mitoticSubRois_1_um'];
    end
else
    subRoiFilename = ip.Results.subRoiFilename;
end

if isempty(ip.Results.subRoiGroupListFilename);
    groupListFilename = subRoiFilename;
else
    groupListFilename = ip.Results.subRoiGroupListFilename;
end

%% Start
groupList1 = groupList;
groupList2= groupList;
groupList1(:,2) = cellfun(@(x) [x filesep 'mtCellEdgeBehaviorClass' filesep subRoiFilename filesep 'sub_1'],groupList(:,2),'uniformoutput',0);
if ip.Results.bipolar
    groupList2(:,2) = cellfun(@(x) [ x filesep 'mtCellEdgeBehaviorClass' filesep subRoiFilename filesep  'sub_2'],groupList(:,2),'uniformoutput',0);
else
    groupList2 = [];
end
groupList = [groupList1;groupList2];

if ip.Results.bipolar && ip.Results.sortGroupList
    btwGrpNames = unique(groupList(:,1),'stable');
    groupListParts = cell(length(btwGrpNames),1,1);
    for i = 1:length(btwGrpNames)
        idxName= find(strcmp(btwGrpNames(i),groupList(:,1)));
        groupListC = groupList(idxName,:);
        % sort by num
        [~,subNum,~] = cellfun(@(x) mitoticUpDirectory(x,5), groupListC(:,2),'uniformoutput',0);
        [~,idxSubNum]  = sort(subNum);
        groupListC = groupListC(idxSubNum,:);
        groupListParts{i}= groupListC;
    end
    
    groupList = vertcat(groupListParts{:});
end

if ~isempty(ip.Results.replaceStr) 
   groupList(:,2) =  cellfun(@(x) strrep(x,ip.Results.replaceStr{1},ip.Results.replaceStr{2}),groupList(:,2),'uniformoutput',0); 
end 

if ~isempty(ip.Results.outputDirectory)
    save([ip.Results.outputDirectory filesep 'groupList_' groupListFilename '.mat'],'groupList'); 
end