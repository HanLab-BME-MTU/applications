function [projList3]=projListSetDiff(projList1,projList2)
% projListSetDiff removes from projList1 the projects in projList2
%
% SYNOPSIS: [projList3]=projListSetDiff(projList1,projList2)
%
% INPUT : projList1/2, two project lists (as from getProj.m)
% OUTPUT: projList3, containing all projects in projList1 that are not in
%         projList2


% format projList1 paths for OS
nProj=length(projList1);
curDir=pwd;
temp1=cellfun(@(x) formatPath(projList1(x,1).anDir),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
temp2=cellfun(@(x) formatPath(projList1(x,1).imDir),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
projList1=cell2struct([temp1 temp2],{'anDir','imDir'},2);

% format projList2 paths for OS
nProj=length(projList2);
curDir=pwd;
temp1=cellfun(@(x) formatPath(projList2(x,1).anDir),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
temp2=cellfun(@(x) formatPath(projList2(x,1).imDir),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
projList2=cell2struct([temp1 temp2],{'anDir','imDir'},2);


temp1=projList2Cell(projList1);
temp1=temp1(:,1);

temp2=projList2Cell(projList2);
temp2=temp2(:,1);

c=1;
rmIdx=[];
for i=1:length(temp2)
    % see if there is a match
    t=find(cellfun(@(x) ~isempty(strmatch(x,temp2{i},'exact')),temp1));
    if ~isempty(t)
        rmIdx=[rmIdx; t];
    end
end

projList3=projList1;
projList3(rmIdx)=[];




