function clickFcn
%get adguiH
adguiH = findall(0,'Tag','adgui');
handles = guidata(adguiH);

clicks = handles.clicks;
handles.clicks = 0;
handles.isFirst = 1;
try
stop(handles.timer)
catch
end

h = warndlg (['in the last 5 seconds you clicked ',num2str(clicks),' times onto the push button'],'Boah!');
set(adguiH,'Visible','off')
uiwait(h)
set(adguiH,'Visible','on')

guidata(adguiH,handles);
