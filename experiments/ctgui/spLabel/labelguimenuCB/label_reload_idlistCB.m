function label_reload_idlistCB
%reloads stored original idlist from labelPanel

%get LP-handle and read old idlist
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
idlist = GetUserData(labelPanelH,'idlist_old');

%overwrite idlist
SetUserData(labelPanelH,idlist,1);

%update labelguiFigure
labelgui('refresh');
