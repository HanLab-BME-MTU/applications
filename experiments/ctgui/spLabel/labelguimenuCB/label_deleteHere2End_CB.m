function label_deleteHere2End_CB(hObject, eventdata, handles)
%function to remove all timepoints from a certain tp on
%based on label_deletetp


cffig = openfig('labelgui','reuse');

imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    return;
end;
dataProperties = GetUserData(imgFigureH,'dataProperties');

ButtonName = questdlg('Remove ALL timepoints from here on?','WARNING','Yes','No','No');
if strcmp(ButtonName,'Yes')
    %continue
else
    return
end

%read data
timeslideH = findall(cffig,'Tag','slider3');
idlist = GetUserData(imgFigureH,'idlist');
curr_time = get(timeslideH,'Value');

%remove all timepoints from here to end
for t = curr_time:length(idlist)
    idlist(t).linklist = [];
end

%write idlist-status
idlist(1).stats.status{end+1}=[date,': deleted frames from ',num2str(curr_time),' to end'];


% %recalculate connections - we don't allow this here
% if recalc
%     idlist = recalcIdlist(idlist,1,[],dataProperties);
% end
        
%save data
SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    view3D_generateHandles;
end

%update labelgui
labelgui('refresh');

