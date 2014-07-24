function label_editlabelistCB
%generates window for editing of labellist (idlist(1).stats.labellist saved in userData of LabelPanel)
%"Calls" label_editlabelCB

edlblH = findall(0,'Type','figure','Name','Label List');
if ~isempty(edlblH)
    delete(edlblH);
end;

idlist = GetUserData(findall(0,'Tag','LabelPanel'),'idlist');
clbllist = idlist(1).stats.labellist;
% add 'new...'
clbllist(end+1) = cellstr('new...');
labellist = char(clbllist(1:end));
cffig = gcbf;
set(cffig,'Units','pixels')
cpos = get(cffig,'Position');

llistfigH = figure('Name','Label List',...
    'Position',[cpos(1)+cpos(3)+5,cpos(2),11*size(labellist,2),15*size(labellist,1)+40],...
    'NumberTitle','off',...
    'MenuBar','None',...
    'Resize','Off');

lboxH = uicontrol('Style','listbox',...
    'Callback','label_editlabelCB',...
    'Position',[6,34,11*size(labellist,2)+10,15*size(labellist,1)],...
    'String',labellist);

lboxH = uicontrol('Style','edit',...
    'Callback','label_editlabelCB',...
    'Position',[6,10,11*size(labellist,2),20],...
    'String',deblank(labellist(1,:)));

