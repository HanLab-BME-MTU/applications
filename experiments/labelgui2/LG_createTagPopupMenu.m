function movieWindowHandles = LG_createTagPopupMenu(movieWindowHandles)
%LG_createTagPopupMenu creates a popup-menu for tag properties

if isempty(movieWindowHandles)
    h = errordlg('movieWindow has been closed','Handle not found');
    uiwait(h);
    return
end

% create a context menu. It has to be stored as a child of the movieWindow,
% which makes it slightly inaccessible. Therefore, we'll store the handle
% also in the guidata structure.
% Since for e.g. lost tags, we can't delete the spot, we need a second
% menu. Create this by copying.
% Actually, it's more flexible to just adjust the menu with the callback
tagPopupMenuH = uicontextmenu('Parent',...
    movieWindowHandles.LG_movieWindow,'Callback','LG_tagPopupMenu_Callback');

% create submenus
% - rename tag (LG_renameTag) 
% - set good tag (LG_setGoodTag)
% - reassign tags (LG_reAssignGUI)
% - delete tag (LG_deleteTag)
% - delete spot (LG_deleteSpot)
% - delete frame (LG_deleteFrame)

% make submenus via callback to tagPopupMenu
renameTagH = uimenu(tagPopupMenuH,'Label','re&name tag');

% set good tag
setGoodTagH = uimenu(tagPopupMenuH,'Label','&set good tag',...
    'Callback','LG_setGoodTag_callback(1)');

% reAssign
uimenu(tagPopupMenuH,'Label','re&link tags',...
    'Callback','LG_reAssignGUI','Separator','on');



% delete tag
uimenu(tagPopupMenuH,'Label','delete &tag',...
    'Callback','LG_deleteTag_callback');

% %--- delete spot is only for some of the menus: Therefore copy menu now
% tagPopupMenuNoSpH = copyobj(tagPopupMenuH,...
%     movieWindowHandles.LG_movieWindow);


% delete spot
uimenu(tagPopupMenuH,'Label','delete &spot',...
    'Callback','LG_deleteSpot_callback');

% delete frame
uimenu(tagPopupMenuH,'Label','delete &frame',...
    'Callback','LG_deleteFrame_callback','Separator','on');

% store handles
movieWindowHandles.tagPopupMenuH = tagPopupMenuH;
%movieWindowHandles.tagPopupMenuNoSpH = tagPopupMenuNoSpH;
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);