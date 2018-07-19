function [groupData,groupValue,groupDataIndividual,groupValueIndividual]=findGroups(list,states,compatible,minLength,minDifferent)
%FINDGROUPS groups consecutive "compatible" entries in the input list
%
%this program is pretty ugly (no reasonable input check etc.). But I needed it real quick...
%
%SYNOPSIS [groupData,groupValue]=findGroups(list,states,compatible)
%
%INPUT    list        : vector to group
%         states      : values in list that should be grouped individually
%         compatible  : (opt) values in states that should be grouped together. 
%                       If several compatible, write e.g. [a,b,c;d,e,NaN]
%         minLength   : (opt) minimum number of entries in a group [1,1] is
%                       length 1 (1 entry)
%         minDifferent: (opt) minimum number of different compatible
%                       entries in a group
%
%OUTPUT   groupData   : [start#, end#] (groups can overlap if compatible-settings are accordingly)
%         groupValue  : values from list (or compatible) corresponding to groups described
%                       by start# and end#
%         *individual : output if only similar values have been grouped
%
%c: jonas, 11/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
if nargin < 3 | isempty(compatible)
    groupMore = 0;
else
    groupMore = 1;
end

if nargin < 4 | isempty(minLength)
    minLength = 1;
end

if nargin < 5 | isempty(minDifferent)
    minDifferent = 1;
end

%if the list is empty: end
if isempty(list)
    groupData = [];
    group = list;
    return
end

nEntries = size(list,1);

%if only one row: end
if nEntries<2
    groupData = [1,1];
    group = list;
    return
end

%end test input

%check groups: group individual if we have to

if nargin > 2 | groupMore == 0
    %calc differences (sum,abs for the matrices)
    deltaPlus = sum(abs(diff(list,1,1)),2);
    deltaMinus = sum(abs(list(1:end-1,:)-list(2:end,:)),2);
    
    startIdx = [1;find(deltaMinus)+1];
    endIdx = [find(deltaPlus);nEntries];
    
    %I assume that startIdx and endIdx must have the same number of entries,
    %but I cannot prove anything
    groupData = [startIdx,endIdx];
    groupValue = list(startIdx,:);
    
    %keep only groups with the right values and the minimum length ([1 1] is length 1)
    goodStatesIdx = find( (ismember(groupValue,states)) & (diff(groupData(:,1:2),1,2)>=minLength-1) );
    
    groupData = groupData(goodStatesIdx,:);
    groupValue = groupValue(goodStatesIdx);
    
    nGroups = size(groupData,1);
    
    %store individual group data
    groupDataIndividual = groupData; 
    groupValueIndividual = groupValue;
end

%group into compatible groups
if groupMore
    %init new groupData
    nComp = size(compatible);
    groupData = [];
    groupValue = [];
    
    %check for occurences of compatible data
    
    for i = 1:nComp(1)
        
        %check that enough values appear in list
        if sum(ismember(compatible(i,:),list))>=minDifferent
            %find where entries in compatible occur
            compList = ismember(list,compatible(i,:));
            %group these entries. We need at least minLength
            [compData] = findGroups(compList,1,[],minLength);
            
            %fill in new data. if there was nothing, we add nothing
            groupData = [groupData;compData];
            if ~isempty(compData)
                %we add the complete compList, even if not everything is
                %present here
                groupValue = [groupValue;repmat(compatible(i,:),[size(compData,1),1])];
            end
        end
    end
    
    %assign output
    [groupData,uniqueIdx] = unique(groupData,'rows');
    groupValue = groupValue(uniqueIdx,:);
    
else
    %do nothing    
end