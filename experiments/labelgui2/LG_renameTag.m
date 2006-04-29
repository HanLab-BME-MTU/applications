function idlist = LG_renameTag(idlist,tagIdx,newTagName)
%LG_renameTag is the function to change a tag label

% change labelcolor
idlist(1).stats.labelcolor{tagIdx} = newTagName;

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: renamed tag %i to %s',date,tagIdx,newTagName);