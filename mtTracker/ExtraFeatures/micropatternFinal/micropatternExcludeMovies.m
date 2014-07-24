function [ micropatternOutput ] = micropatternExcludeMovies(micropatternOutput,maskID,saveDir)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%

if isempty(micropatternOutput) || nargin <1
    [file, path]   = uigetfile(pwd,'Please Select micropatternOutput.mat file');
    s = load([path filesep file]);
    micropatternOutput = s.micropatternOutput;
end

if isempty(maskID) || nargin <2 % get rid of function for all masks;
    maskID = 1:numel(micropatternOutput);
    maskID = maskID';
end

if isempty(saveDir) || nargin<3
    saveDir = uigetdir(pwd,'Where would you like to save the modified micropatternOutput?');
end

%%
% loop through all masks



% load groupListWholeCell
groupList = micropatternOutput{1}.groupListWholeCell ;


choices = groupList(:,2);

% let the user pick a cell to remove
userPick = menu('please select a cell to remove',choices);


toRemove = choices(userPick); % get the name to remove
groupList(userPick,:) = [];
groupListWholeCell = groupList; 
%micropatternOutput.groupListWholeCell = groupList;



% loop through all the masks find the cell and remove
for iMask = 1:length(maskID)
    
    if isfield(micropatternOutput{iMask},'removedByUser');
        numPrevRemoved = numel(micropatternOutput{maskID(iMask)}.removedByUser);
        micropatternOutput{maskID(iMask)}.removedByUser{numPrevRemoved+1,1} = toRemove;
    else
        micropatternOutput{iMask}.removedByUser{1} = toRemove;
    end
    clear groupList
    clear idx
    for iExtract = 1:numel(micropatternOutput{maskID(iMask)}.maskParams.extractType)
        
        if micropatternOutput{maskID(iMask)}.maskParams.subRegions ==1
            
            for iReg = 1:numel(micropatternOutput{maskID(iMask)}.groupList{1})
                for iWind = 1:numel(micropatternOutput{maskID(iMask)}.groupList{1}{1})
                    groupList = micropatternOutput{maskID(iMask)}.groupList{iExtract}{iReg}{iWind} ;% maybe change fieldname
                    
                    idx = 1:length(groupList(:,2));
                    % Remove the user selected files ( i know there is an
                    % easier way to do this :(
                    for iProj = 1:length(groupList(:,2))
                        
                        if strfind(char(groupList(iProj,2)),char(toRemove)) == 1;
                            idx(iProj) = 1;
                        else
                            idx(iProj) = 0 ;
                        end
                    end
                    groupList(logical(idx'),:) = [];
                    % save the groupList
                    micropatternOutput{maskID(iMask)}.groupList{iExtract}{iReg}{iWind} = groupList;
                    micropatternOutput{maskID(iMask)}.groupListWholeCell = groupListWholeCell;
                    
                end
            end
            
        else  % a bit redundant if I was smarted with my indexing could maybe have
            % saved some code.
            
            for iWind = 1:numel(micropatternOutput{maskID(iMask)}.groupList{1})
                groupList = micropatternOutput{maskID(iMask)}.groupList{iExtract}{iWind};
                
                % Remove the user selected files
                for iProj = 1:length(groupList(:,2))
                    
                    if strfind(char(groupList(iProj,2)),char(toRemove)) == 1
                        
                        idx(iProj) = 1;
                    else
                        idx(iProj) = 0 ;
                    end
                    
                end % for iproj
                
                groupList(logical(idx'),:) = [];
                
                
                
                
                
                
                
                % save the new groupList
                micropatternOutput{maskID(iMask)}.groupList{iExtract}{iWind} = groupList;
                
            end
            
        end
        
    end % for iExtract
end % for iMask

save([saveDir filesep 'micropatternOutput.mat'], 'micropatternOutput');

end

