function [ groupList] = makeCorticalGroupList( groupList ,name,bipolar)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
projList = groupList(:,2); 
if bipolar ==1
 for iProj = 1:length(groupList(:,2))
 
 anDir = char(projList(iProj)); 
s  = load([anDir filesep 'poleInfo.mat']);
%     if isfield(s.poleInfo,'mult') % this was to fix the fact that I used a slightly different format for poleInfo in 
%            % the 2-11 data. dipshit... 
%            s.poleInfo.numPoles = 2; 
%        end
    
    if s.poleInfo.numPoles ==2
        projListFlag(iProj,1)=1;
    else
        projListFlag(iProj,1)=0;
    end

 end 
else projListFlag = ones(length(groupList(:,2)),1) ; 
    
end
groupList = groupList(logical(projListFlag),:); % only take bipolar
    groupList1 = groupList; 
    groupList2= groupList; 
 groupList1(:,2) = cellfun(@(x) strrep(x,'roi_1',['roi_1' filesep name filesep 'sub_1']),groupList(:,2),'uniformoutput',0); 
 if bipolar == 1
 groupList2(:,2) = cellfun(@(x) strrep(x,'roi_1',['roi_1' filesep name filesep  'sub_2']),groupList(:,2),'uniformoutput',0); 
 else 
     groupList2 = []; 
 end 
 groupList = [groupList1;groupList2];
% sort 
%names = unique(groupList(:,1),'stable'); 

% for i = 1:length(names)
% idx{i} = find(strcmp(names{i},groupList(:,1))); 
% end

%save('subRois_Cortical_All.mat','groupList'); 
end