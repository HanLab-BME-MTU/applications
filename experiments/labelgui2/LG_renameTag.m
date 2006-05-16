function idlist = LG_renameTag(idlist,tagIdx,newTagName)
%LG_renameTag is the function to change a tag label

% check for -new-
if strcmp(newTagName,'-new-')
    newTagName = inputdlg('Please enter new name');
    if isempty(newTagName) || strcmp(newTagName{1},'-new-')
        % don't assign anything if user cancelled
    else
        newTagName = newTagName{1};
        % add to labellist
        idlist(1).stats.labellist{end+1} = newTagName;
    end
end

% change labelcolor
idlist(1).stats.labelcolor{tagIdx} = newTagName;

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: renamed tag %i to %s',date,tagIdx,newTagName);