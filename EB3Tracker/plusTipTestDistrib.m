function [discrimMat]=plusTipTestDistrib(groupData,distribs2test)
% plusTipTestDistib returns discrimination matrices for plus tip distributions
%
% SYNOPSIS : [discrimMat]=plusTipTestDistrib(groupData,distribs2test)
%
% INPUT : 
%   groupData :     nGroups x 1 structure containing fields for
%   distributions
%                   from plus tip data - gs,fs,bs,gl,...bd (see
%                   plusTipPoolGroupData)
%   distribs2test : string or cell array of multiple strings containing
%                   distribution identifier(s):
%                       'gs' growth speed 'fs' fgap speed 'bs' bgap speed
%                       'gl' growth lifetime 'fl' fgap lifetime 'bl' bgap
%                       lifetime 'gd' growth displacement 'fd' fgap
%                       displacement 'bd' bgap displacement
%
% OUTPUT :
%   disrimMat : structure containg matrices with p-values above the
%               diagonal from the permutation test for the means, and with
%               percent confidence below the diagonal from the bootstrapped
%               distribution test, which gives the percent confidence that
%               the two distributions are different.
%
%               the fields with the suffix "cell" are identical to the
%               matrix form except they also include group labels in case
%               you want to write them out into an Excel file.
%  

nGroups=length(groupData);

% keep track of the group names in order
temp=struct2cell(groupData);
grpNames=cell(nGroups,1);
for i=1:nGroups
    grpNames{i,1}=temp{1,i}.name;
end

% remove .info field and get fieldnames corresponding to distributions
groupData=rmfield(groupData,'info');
fN=fieldnames(groupData);


% figure out which distributions the user wants to test
if nargin<2 || isempty(distribs2test)
    idx=(1:9)'; % do them all
else
    if ischar(distribs2test)
        distribs2test={distribs2test};
    end
    if iscell(distribs2test)
        idx=zeros(length(distribs2test),1);
        for iDist=1:length(distribs2test)
            idx(iDist)=find(cellfun(@(x) ~isempty(strmatch(distribs2test{iDist},x,'exact')),fN));
        end
    else
        error('distribs2test must be a string or a cell');
    end
end




for i=1:length(idx)
    for iGroup=1:nGroups
        dataStruct(iGroup,1).(fN{idx(i)})=groupData(iGroup,1).(fN{idx(i)});
    end
    % do perm test below and distrib test above the diag
    testStruct.(fN{idx(i)})=[21 20];
end


% run the test
discrimMat=discriminationMatrix(dataStruct,testStruct);

% convert the matrices to cell arrays with labels
for i=1:length(idx)

    switch fN{idx(i)}
        case 'gs'
            str='growth speed';
        case 'fs'
            str='fgap speed';
        case 'bs'
            str='bgap speed';
        case 'gl'
            str='growth lifetime';
        case 'fl'
            str='fgap lifetime';
        case 'bl'
            str='bgap lifetime';
        case 'gd'
            str='growth displacement';
        case 'fd'
            str='fgap displacement';
        case 'bd'
            str='bgap displacement';
    end

    % convert from matrix to cell array so we can pad with labels
    temp=mat2cell(discrimMat.(fN{idx(i)}),ones(nGroups,1),ones(nGroups,1));
    % add stat name in upper left and labels along rows and columns
    discrimMat.([fN{idx(i)} '_cell'])=[[{str} grpNames'];[grpNames temp]];

end


