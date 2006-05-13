function LG_tagPopupMenu_Callback
% this is the callback that is executed just before the tagPopupMenu is launched

% Here, we add the submenu with the tagNames and enable/disable menu items
% if necessary
% gco returns the object that has been clicked on
% gcbo returns the handle to the popupMenu

% careful! Changes we make to the menu will persist!
% menuOrder (reversed)
% - rename tag (LG_renameTag) -- end
% - set good tag (LG_setGoodTag)
% - reassign tags (LG_reAssignGUI)
% - delete tag (LG_deleteTag)
% - delete spot (LG_deleteSpot)
% - delete frame (LG_deleteFrame) -- 1

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;
tagPopupMenuH = gcbo;
tagPopupMenuItemH = get(tagPopupMenuH,'Children');
objectH = gco;
tagIdx = get(objectH,'UserData');
tagLabel = movieWindowHandles.idlistData.labelcolor{tagIdx};
currentTime = LG_getCurrentTime;

% get handle to the rename-menuItem. It's the first item, so it's last in
% the list of children
renameTagH = tagPopupMenuItemH(end);
labellist = movieWindowHandles.idlistData.labellist;
linklist = movieWindowHandles.idlist(currentTime).linklist;

% delete previous children first!
delete(get(renameTagH, 'Children'))

for i=1:size(labellist,1)
    uih=uimenu(renameTagH,'Label',labellist{i},...
        'Callback',['LG_renameTag_callback(''',labellist{i},''')']);
    if strcmp(labellist{i},tagLabel)
        set(uih,'Checked','on')
    end
end

% check whether the tag is already a 'good' tag. If yes, allow to set bad
% tag
if linklist(tagIdx,5) == 2
    setGoodTagH = tagPopupMenuItemH(end-1);
    set(setGoodTagH, 'Label','&set good tag',...
        'Callback','LG_setGoodTag_callback(true)');
else
    setGoodTagH = tagPopupMenuItemH(end-1);
    set(setGoodTagH, 'Label','&set bad tag',...
        'Callback','LG_setGoodTag_callback(false)');
end


% check whether (a) the current frame is a fusion, and (b) whether the
% parent of the current object is axes or figure.
% From the info-figure, we can delete a tag out of a fusion, from the
% movieWindow, we can delete the entire spot, but not the other way round.
% Also, we can't rename a fusion.
parentH = get(objectH,'Parent');
parentIsAxes = strcmp(get(parentH,'Type'),'axes');

isFusion = any(linklist(tagIdx,3) == [3,4,5]);
deleteSpotH = tagPopupMenuItemH(end-4);
deleteTagH = tagPopupMenuItemH(end-3);
% we already know renameTagH

if parentIsAxes
    if isFusion
        % hide rename, deleteTag
        set(deleteSpotH,'Visible','on');
        set([deleteTagH; renameTagH; setGoodTagH],'Visible','off');
    else
        % make all visible
        set([deleteSpotH;deleteTagH; renameTagH; setGoodTagH],'Visible','on');
    end
else % parent is window

    % hide deleteTag
    set(deleteTagH,'Visible','off');
    set([deleteSpotH; renameTagH; setGoodTagH],'Visible','on');
end



% check whether it's an estimated spot - no deleteSpot!
if linklist(tagIdx,3) == 1
    set(deleteSpotH,'Visible','off');
else
    set(deleteSpotH,'Visible','on');
end


