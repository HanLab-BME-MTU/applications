function view3D_refresh
%updates  figure: show spots and arrows


%-------------get necessary handles and data
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
view3DH = GetUserData(labelPanelH,'view3DGUIH');
if isempty(view3DH)
    error('this is a bug which still needs fixing')
end
view3D_handles = guidata(view3DH); %handle structure of view3D
idlist = GetUserData(labelPanelH,'idlist');
labelguiH = findall(0,'Tag','labelgui');
label_handles = guidata(labelguiH); %handle structure of labelGui
view3DfigH = GetUserData(labelPanelH,'view3DfigH');

%plot data
lastLabelguiTime = view3D_handles.lastLabelguiTime;
arrowHandles = view3D_handles.arrowHandles;
colorList = view3D_handles.colorList;
nspots = view3D_handles.nspots;
%view3Dgui data
drawLinks = get(view3D_handles.view3D_link_PD,'Value');
timeSetBy = get(view3D_handles.view3D_setTime_PD,'Value');
tstart = get(view3D_handles.view3D_tstart_slider,'Value');
tend = get(view3D_handles.view3D_tend_slider,'Value');
%labelgui data
labelguiTime = get(label_handles.slider3,'Value');
% %insert here: labelgui stack data

%-----------adjust time
switch timeSetBy
    case 1 %independent
        %turn off last mark
        lastLabelguiTimeDataHandles = arrowHandles{lastLabelguiTime};
        for i = 1:nspots(lastLabelguiTime)
            set(lastLabelguiTimeDataHandles(i,3),'EdgeColor','none');
        end
    case 2 %time set by labelgui
        if labelguiTime>tend
            set(view3D_handles.view3D_tend_slider,'Value',labelguiTime);
            set(view3D_handles.view3D_tend_txt,'String',num2str(labelguiTime));
            tend = labelguiTime;
        elseif labelguiTime<tstart
            set(view3D2_handles.tstart_slider,'Value',labelguiTime);
            set(view3D2_handles.tstart_txt,'String',num2str(labelguiTime));
            tstart = labelguiTime;
        end
        labelguiTime = get(label_handles.slider3,'Value');
        %mark labelguiTime, turn off last mark
        lastLabelguiTimeDataHandles = arrowHandles{lastLabelguiTime};
        for i = 1:nspots(lastLabelguiTime) %turn off
            set(lastLabelguiTimeDataHandles(i,3),'EdgeColor','none');
        end
        if ~isempty(idlist(labelguiTime).linklist) %set
            labelguiTimeDataHandles = view3D_handles.arrowHandles{labelguiTime};
            set(labelguiTimeDataHandles(:,3),'EdgeColor','k');
            view3D_handles.lastLabelguiTime = labelguiTime;%only update lastLabelguiTime if new marker has been set
        end
    case 3 %tend->labelguiTime
        set(label_handles.slider3,'Value',tend);
        labelgui('slider3_Callback',label_handles.slider3,[],label_handles);
        labelguiTime = get(label_handles.slider3,'Value');
        %mark labelguiTime, turn off last mark
        lastLabelguiTimeDataHandles = arrowHandles{lastLabelguiTime};
        for i = 1:nspots(lastLabelguiTime) %turn off
            set(lastLabelguiTimeDataHandles(i,3),'EdgeColor','none');
        end
        if ~isempty(idlist(labelguiTime).linklist) %set
            labelguiTimeDataHandles = view3D_handles.arrowHandles{labelguiTime};
            set(labelguiTimeDataHandles(:,3),'EdgeColor','k');
            view3D_handles.lastLabelguiTime = labelguiTime; %only update lastLabelguiTime if new marker has been set
        end
end



%-----------show data

%turn everything off
set(cat(1,arrowHandles{1:end}),'Visible','off');

%turn on spots and arrows and arrowheads tstart:tend
arrowHandles2Draw = cat(1,arrowHandles{tstart:tend});
colorList2Draw = cat(1,colorList{tstart:tend});
switch (drawLinks-1)*(tend>tstart)+1 % = 1 if tend==tstart, else show what user selected
    case 1 %none at all
        if ~isempty(arrowHandles2Draw)
            set(arrowHandles2Draw(:,3),'Visible','on');
        end
    case 2 %black lines
        if ~isempty(arrowHandles2Draw)
            set(arrowHandles2Draw(:,3),'Visible','on');
            set(arrowHandles2Draw(:,1),'Visible','on','Color','k');
        end
    case 3 %colored lines
        if ~isempty(arrowHandles2Draw)
            set(arrowHandles2Draw(:,3),'Visible','on');
            for i = 1:size(colorList2Draw,1)
                set(arrowHandles2Draw(i,1),'Visible','on','Color',colorList2Draw(i,:));
            end
        end
    case 4 %black arrows
        if ~isempty(arrowHandles2Draw)
            set(arrowHandles2Draw(:,3),'Visible','on');
            set(arrowHandles2Draw(:,1),'Visible','on','Color','k');
            set(arrowHandles2Draw(:,2),'Visible','on','FaceColor','k');
        end
    case 5 %colored arrows
        if ~isempty(arrowHandles2Draw)
            set(arrowHandles2Draw(:,3),'Visible','on');
            for i = 1:size(colorList2Draw,1)
                set(arrowHandles2Draw(i,1),'Visible','on','Color',colorList2Draw(i,:));
                set(arrowHandles2Draw(i,2),'Visible','on','FaceColor',colorList2Draw(i,:));
            end
        end
    case 6 %colored lines only
        if ~isempty(arrowHandles2Draw)
            for i = 1:size(colorList2Draw,1)
                set(arrowHandles2Draw(i,1),'Visible','on','Color',colorList2Draw(i,:));
            end
        end
    case 7 %colored arrows only
        if ~isempty(arrowHandles2Draw)
            for i = 1:size(colorList2Draw,1)
                set(arrowHandles2Draw(i,1),'Visible','on','Color',colorList2Draw(i,:));
                set(arrowHandles2Draw(i,2),'Visible','on','FaceColor',colorList2Draw(i,:));
            end
        end
end

%plot axis
if get(view3D_handles.view3D_check_drawAxis,'Value')==1
    axisCoordinates = view3D_handles.axisCoordinates;
    figure(view3DfigH);
    delete(view3D_handles.lineH);
    lineH = line(squeeze(axisCoordinates(:,1,tstart:tend)),squeeze(axisCoordinates(:,2,tstart:tend)),squeeze(axisCoordinates(:,3,tstart:tend)),'Color','k');
    view3D_handles.lineH = lineH; %overwrite existing lineHandle
elseif ~isempty(view3D_handles.lineH)
    set(view3D_handles.lineH,'Visible','off');
end

%save data
guidata(view3D_handles.view3DGUI,view3D_handles);


%make gui and figure topmost
figure(view3DH);
figure(view3DfigH);
daspect([1,1,1]);