function manrelink
% manrelink creates objects with callbacks to manipulate cell tracks
%
% SYNOPSIS       manrelink
%
% INPUT          none (gets data from the current object, which is the object with the callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   manrelink uses {nothing}
%                                  
%                manrelink is used by { changeframe (as a callback) }
%
% Colin Glass, Feb 04        

%this is the callback of the cell numbers in the shown image(look at
%changeframe). 

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

%first we have to find out which cell the user has clicked on
whichcell=get(gco,'String')


% posi= get(gco,'Position');
% where1=posi(1);
% where2=posi(2);
% 
% pospic=get(gcf,'Position');
% wherepic1=pospic(1);
% wherepic2=pospic(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program creates three objects:
%- listbox(linklist) where the program keeps track of the cells you have selected for
%  linkage. When you click on one of the listed cells, it will be linked to
%  the current callback cell
%- pushbutton (linkbutton). If you click on this, the current callback cell
%  will be added to the list (in listbox)
%- pushbutton (chuckout). If you click on this, the current callback cell
%  will be erased
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first we look if these three bjects already exist somewhere and, if that
%should be the case, erase them.
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



%get the current list that should go into listbox
listofcells=handles.listofcells;

%find the position of the users click, so that we can create the objects
%there
figure(gcbf);
%position whithin figure
pos = get(gcf,'currentPoint')
%position of figure
imgPos = get(gca,'Position');

%calculate the position the objects should have
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


%create the objects with the respective callbacks

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


   
   
   