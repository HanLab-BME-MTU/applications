function manrelink
hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

whichcell=get(gco,'String')

posi= get(gco,'Position');
where1=posi(1);
where2=posi(2);

pospic=get(gcf,'Position');
wherepic1=pospic(1);
wherepic2=pospic(2);


linklistH =findall(0,'Style','listbox','Tag','linklist');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');
chuckoutH =findall(0,'Style','pushbutton','Tag','chuckout');

if ~isempty(linklistH)
    delete(linklistH);
end;

if ~isempty(linkbuttonH)
    delete(linkbuttonH);
end;
if ~isempty(chuckoutH)
    delete(chuckoutH);
end;

listofcells=handles.listofcells;
figure(gcbf);
pos = get(gcf,'currentPoint')
% if strcmp(get(gca,'Units'),'normalized')
%     set(gca,'Units','pixels');
% end
imgPos = get(gca,'Position');
%adjust position of poplistbox so that its top is not outside current figure
lbPos = [pos(1)  pos(2) 60 60]
linkaddbuttonpos=[pos(1)  pos(2)+60 60 30]
delbuttonPos = [pos(1)  pos(2)+90 60 30]

% % deltaH = (buttonPos(2)+buttonPos(4))-(imgPos(4)); %height of box+button minus imgHeight in figure coords
% % if deltaH>0 %if position of button outside img
% %     %adjust poplistbox/button-position
% %     plbPos(2) = plbPos(2) - deltaH;
% %     buttonPos(2) = buttonPos(2) - deltaH;
% % end
% % %make the same for the side. This time, check with figure size
% % figPos = get(gcbf,'Position');
% % deltaR = plbPos(1)+plbPos(3)-figPos(3);
% % if deltaR>0
% %     plbPos(1) = plbPos(1) - deltaR;
% %     buttonPos(1) = buttonPos(1) - deltaR;
% % end
% % 


linklistH = uicontrol('Style', 'listbox',...
    'Callback','linknow',...
    'String', listofcells,...
    'UserData',whichcell,...
    'Tag','linklist',...
    'Position',lbPos)%,...
    %'Position',[300,100,100,100]);


linkbuttonH = uicontrol('Style', 'pushbutton', 'String', 'Link',...
    'Callback','linkhelp',...
    'Tag','linkbutton',...
    'UserData',whichcell,...
    'Position',linkaddbuttonpos)%,...
    %'Position', [200,100,100,100]);


chuckoutH = uicontrol('Style', 'pushbutton', 'String', 'Delete',...
    'Callback','deletingrow',...
    'Tag','chuckout',...
    'UserData',whichcell,...
    'Position',delbuttonPos)%,...
   %'Position', [100,100,100,100]);


   
   
   