function LG_navi_menuLoadAll(hObject, eventdata, handles)
%LG_navi_menuLoadAll is the main loader function for labelgui2

% load movie, launch movie window
success = LG_navi_menuLoadMovie(0,0,handles);

% if there was no success: die (we warned before)
if ~success
    return
end

% load idlist. Reload handles first.
handles = LG_getNaviHandles;
success = LG_navi_menuLoadIdlist(0,0,handles);

% if no success: We need to call the plotting function from here
if ~success
    LG_gotoFrame(1);
end

