function label_showIntFig(hObject,eventdata,handles)
%if exist, loads and shows intFigure for current idlist

%is item checked
isChecked = get(hObject,'Checked');
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');


if strcmp(isChecked,'on') %uncheck, close all windows
	set(hObject,'Checked','off');
	intFigH = findall(0,'Tag','intFig');
	close(intFigH);
    %delete all the saved intFigH
    labelPanelList = findall(0,'Tag','LabelPanel');
    for i = 1:length(labelPanelList)
        SetUserData(labelPanelList(i),[],1,'intFigH');
    end
	return
	
else %check and open figure
	set(hObject,'Checked','on');
end

%get data
labelguiH = handles.labelgui;

%where is the figure stored?
dataFile = GetUserData(imgFigureH,'dataFile');

if isempty(dataFile)
    [fname,pathName] = uigetfile({'intFigure??-???-????-??-??-??.fig','intFigures'},'select intFigure file');
    if(fname(1)==0)
        image = [];
        fname = [];
        return;
    end;
    
    
else
    
    pathName = dataFile.path;
    fname=chooseFile('intFigure',pathName,'new');
    if isempty(fname)
        [fname,pathName] = uigetfile({'intFigure??-???-????-??-??-??.fig','intFigures'},'select intFigure file');
        if(fname(1)==0)
            image = [];
            fname = [];
            return;
        end;
    end
end

cd(pathName);
figH = openfig(fname,'new'); 
set(figH,'Tag','intFig');

%if the line handle to the intensity lines has been stored, display a
%button that allows to update this
figHandles = guidata(figH);
figure(figH);
lineHCandidates = findall(figH,'Type','line');
if isempty(figHandles)
    %if the handles have not been stored, try to find the correct lines
    %other lines are blue,red,magenta and black
    lineColors = cell2mat(get(lineHCandidates,'Color'));
    [dummy,goodColorIdx] = setdiff(lineColors,[0 0 0;0 0 1;1 0 0;1 0 1],'rows');
    figHandles.spotIntLineH = lineHCandidates(goodColorIdx);
    guidata(figH,figHandles);
end

if ~isempty(figHandles)
    uh = uicontrol('Style','pushbutton',...
        'Tag','updateIntFig_PB',...
        'Position',[5,5,120,23],...
        'Callback','label_updateIntFigCB(gcbo,[],guidata(gcbo))','String','update tag intensities');
    set(uh,'Units','characters');
end

%set buttonDownFcn, so that a click on the figure changes the current
%timepoint of labelgui
set(gca,'ButtonDownFcn','label_gotoFrame_BDFCN');
%do it also for all lines
set(lineHCandidates,'ButtonDownFcn','label_gotoFrame_BDFCN');

if ishandle(imgFigureH)
    figure(imgFigureH);
end
figure(labelguiH);

%save handle
SetUserData(imgFigureH,figH,1,'intFigH');