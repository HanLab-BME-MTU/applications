function [discrimMat]=plusTipTestDistrib(saveDir,groupData,distribs2test)
% plusTipTestDistib returns discrimination matrices for plus tip distributions
%
% SYNOPSIS : [discrimMat]=plusTipTestDistrib(groupData,distribs2test)
%
% INPUT : 
%   groupData :     nGroups x 1 structure containing fields for
%                   distributions from plus tip data - gs,fs,bs,gl,...bd 
%                   (see plusTipPoolGroupData)
%   distribs2test : string or cell array of multiple strings containing
%                   distribution identifier(s):
%                       'gs' growth speed 'fs' fgap speed 'bs' bgap speed
%                       'gl' growth lifetime 'fl' fgap lifetime 'bl' bgap
%                       lifetime 'gd' growth displacement 'fd' fgap
%                       displacement 'bd' bgap displacement
%
% OUTPUT :
%   disrimMat : structure containing matrices with p-values above the
%               diagonal from the permutation test for the means, and with
%               percent confidence below the diagonal from the bootstrapped
%               distribution test, which gives the percent confidence that
%               the two distributions are different.
%
%               the fields with the suffix "cell" are identical to the
%               matrix form except they also include group labels in case
%               you want to write them out into an Excel file.
%  

if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select Output Directory.');
end



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

 
 hitList=0;
 
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
    stringency = 0.05;

    hitsIdx = find(discrimMat.(fN{idx(i)})(1,:) <stringency);
    
   
    
    
    if hitList == 0;
    for iHit = 1:length(hitsIdx)
        hits{iHit,1} = grpNames((hitsIdx(iHit)),1);
        hits{iHit,2} = fN{idx(i)};
        hitList = 1;
    end 
        else
    
        for iHit =1:length(hitsIdx)
            curLength = length(hits(:,1));
            hits{curLength+1,1} = grpNames((hitsIdx(iHit)),1);
            hits{curLength+1,2} = fN{idx(i)};
        end
   
    end  
     matLength = length(discrimMat.([fN{idx(i)} '_cell'])(1,:));
     
     space = cell(1,matLength);
    if i ==1 
    discrimMat.All = [discrimMat.([fN{idx(i)} '_cell']); space]; 
    else 
        discrimMat.All =[discrimMat.All ;discrimMat.([fN{idx(i)} '_cell'])];
        discrimMat.All = [discrimMat.All; space];
    end 
   

end

if hitList == 0
    hits = 'No Significant Differences';
else 
end 
discrimMat.hits = hits;
fileName = ['discrimMat','.mat'];
save([saveDir filesep fileName],'discrimMat');







end  

