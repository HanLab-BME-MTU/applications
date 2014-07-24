function popuphere
%generates popup window to change spot label. "Calls" poplistboxCB

idlist=GetUserData(GetUserData(openfig('labelgui','reuse'),'currentWindow'),'idlist');
clabellist=(idlist(1).stats.labellist);
labellist=char(clabellist(1:end));

hpop=findall(0,'Style','listbox','Tag','selcetlabel');
hpush=findall(0,'Style','pushbutton','Tag','deletespot');

if ~isempty(hpop)
    delete(hpop);
end;
if ~isempty(hpush)
    delete(hpush);
end;
figure(gcbf);

pos = get(gcf,'currentPoint');
if strcmp(get(gca,'Units'),'normalized')
    set(gca,'Units','pixels');
end
imgPos = get(gca,'Position');
%adjust position of poplistbox so that its top is not outside current figure
plbPos = [pos(1)  pos(2) 15*size(labellist,2) 18*size(labellist,1)];
buttonPos = [pos(1)  pos(2)+18*size(labellist,1) 15*size(labellist,2) 20];
deltaH = (buttonPos(2)+buttonPos(4))-(imgPos(4)); %height of box+button minus imgHeight in figure coords
if deltaH>0 %if position of button outside img
    %adjust poplistbox/button-position
    plbPos(2) = plbPos(2) - deltaH;
    buttonPos(2) = buttonPos(2) - deltaH;
end
%make the same for the side. This time, check with figure size
figPos = get(gcbf,'Position');
deltaR = plbPos(1)+plbPos(3)-figPos(3);
if deltaR>0
    plbPos(1) = plbPos(1) - deltaR;
    buttonPos(1) = buttonPos(1) - deltaR;
end


hpop = uicontrol('Style', 'listbox',...
    'Callback','poplistboxCB',...
    'String', labellist,...
    'UserData',gcbo,...
    'Tag','selcetlabel',...
    'Position', plbPos);

hpush = uicontrol('Style', 'pushbutton', 'String', 'Delete',...
    'Callback',['label_deleteSpotCB(' num2str(get(gcbo,'UserData')) ')'],...
    'Tag','deletespot',...
    'UserData',gcbo,...
    'Position', buttonPos);

sel = strmatch(get(gcbo,'String'),get(hpop,'String'),'exact');
if~isempty(sel);
    set(hpop,'Value',sel);
end;