function windowHandles = LG_createTagPopupMenu(windowHandles,figureHandle)
%LG_createTagPopupMenu creates a popup-menu for tag properties

if ~ishandle(figureHandle)
    h = errordlg('movieWindow has been closed','Handle not found');
    uiwait(h);
    return
end

% create a context menu. It has to be stored as a child of the respective,
% which makes it slightly inaccessible. 
% Since for e.g. lost tags, we can't delete the spot, we need a second
% menu. Create this by copying.
% Actually, it's more flexible to just adjust the menu with the callback
tagPopupMenuH = uicontextmenu('Parent',...
    figureHandle,'Callback','LG_tagPopupMenu_Callback',...
    'SelectionHighlight','on');

% create submenus
% - rename tag (LG_renameTag)
% - change color (LG_setTagColor)
% - set good tag (LG_setGoodTag)
% - reassign tags (LG_reAssignGUI)
% - delete tag (LG_deleteTag)
% - delete spot (LG_deleteSpot)
% - delete frame (LG_deleteFrame)

% make submenus via callback to tagPopupMenu
uimenu(tagPopupMenuH,'Label','re&name tag');

% change color
uimenu(tagPopupMenuH,'Label','change &color',...
    'Callback','LG_changeTagColor');

% set good tag
uimenu(tagPopupMenuH,'Label','&set good tag',...
    'Callback','LG_setGoodTag_callback(1)','Separator','on');

% reAssign
uimenu(tagPopupMenuH,'Label','re&link tags',...
    'Callback','LG_reAssignGUI','Separator','on');



% delete tag
uimenu(tagPopupMenuH,'Label','delete &tag',...
    'Callback','LG_deleteTag_callback');


% delete spot
uimenu(tagPopupMenuH,'Label','delete &spot',...
    'Callback','LG_deleteSpot_callback');

% delete frame
uimenu(tagPopupMenuH,'Label','delete &frame',...
    'Callback','LG_deleteFrame_callback','Separator','on');

% store handles
windowHandles.tagPopupMenuH = tagPopupMenuH;