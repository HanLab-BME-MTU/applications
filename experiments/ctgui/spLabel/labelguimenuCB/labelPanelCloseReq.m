function labelPanelCloseReq
%function for gracefully closing the labelpanel

%find labelPanel: copied from matlabs closereq
shh=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
imgFigureH=get(0,'CurrentFigure');
set(0,'ShowHiddenHandles',shh);

try
    
    %delete all associated windows
    
    %3d-viewer
    view3DH = GetUserData(imgFigureH,'view3DGUIH');
    if ishandle(view3DH)
        %remember position
        positions.view3DPos = get(view3DH,'Position');
        delete(view3DH);
    end
    view3DfigH = GetUserData(imgFigureH,'view3Dfig');
    if ishandle(view3DfigH)
        %remember position
        positions.view3DfigPos = get(view3DfigH,'Position');
        delete(view3DfigH);
    end
    
    %xz/yz-viewer
    zFigH = GetUserData(imgFigureH,'XZYZFigureH');
    if ishandle(zFigH)
        %remember position
        positions.zFigPos = get(zFigH,'Position');
        delete(zFigH);
    end
    
    %movieDataFigure
    dataFigH = GetUserData(imgFigureH,'CurrentMovieData');   
    if ishandle(dataFigH)
        %remember position
        positions.dataFigPos = get(dataFigH,'Position');
        delete(dataFigH);
    end
    
    %intFigure
    figH = GetUserData(imgFigureH,'intFigH');
    if ishandle(figH)
        %remember position
        positions.intFigPos = get(figH,'Position');
        delete(figH);
    end
    
    %disFigure
    figH = GetUserData(imgFigureH,'disFigH');
    if ishandle(figH)
        %remember position
        positions.disFigPos = get(figH,'Position');
        delete(figH);
    end
    
    %relinkgui
    relinkH = findall(0,'Tag','reLinkGUI');
    if ishandle(relinkH)
        positions.relinkPos = get(relinkH,'Position');
        delete(relinkH);
    end
    
    %remember labelPanelPosition
    positions.labelPanelPos = get(imgFigureH,'Position');
    
    %remove entry from labelgui-pd
    labelguiH = findall(0,'Tag','labelgui');
    
    %if the labelgui has already been closed: close without making a fuss
    if isempty(labelguiH)
        %use delete, otherwise, the function calls itself recursively
        delete(imgFigureH)
        return
    end
    
    %save positions
    positionNames = fieldnames(positions);
    positionValues = struct2cell(positions);
    %retrieve previously saved pos
    savedPositions = GetUserData(labelguiH,'positions');
    for nPos = 1:length(positionNames)
        eval(['savedPositions.',positionNames{nPos},'=',mat2str(positionValues{nPos}),';']);
    end
    %save in labelgui
    SetUserData(labelguiH,savedPositions,1,'positions');
    
    labelguiHandles = guidata(labelguiH);
    
    pdH = labelguiHandles.label_chooseWindow_PD;
    
    %to find out which entry to remove: get handle list
    handleList = GetUserData(labelguiH,'labelPanelHandles');
    
    myHandleIdx = find(handleList == imgFigureH);
    
    %get pd-string
    pdString = get(pdH,'String');
    
    %delete entry from pdString & from handleList
    pdString(myHandleIdx) = [];
    handleList(myHandleIdx) = [];
    
    
    %check whether labelpanel being deleted is the active labelPanel
    activeLabelPanelH = GetUserData(labelguiH,'currentWindow');
    if activeLabelPanelH == imgFigureH
        %we are deleting the active figure! Set pd to the first entry or write
        %'no movie loaded'
        if isempty(handleList)
            pdString = [{'no movie loaded'};pdString];
            currentWindow = [];
        else
            currentWindow = handleList(1);
        end
        %set correct currentWindow
        SetUserData(labelguiH,currentWindow,1);
        set(pdH,'Value',1);
        
        %we can call refresh already here, because it only looks at
        %currentwindow
        labelgui('refresh');
        
    else %make sure that we choose the right menu from the list
        activeIdx=find(activeLabelPanelH == handleList);
        set(pdH,'Value',activeIdx);
    end
        
       
    
    
    %save entries
    set(pdH,'String',pdString);
    SetUserData(labelguiH,handleList,1,'labelPanelHandles');
    
catch
    errmsg = lasterr;
    h = warndlg(['error: ',errmsg,' window will close anyway']);
    uiwait(h)
end

%kill labelPanel
delete(imgFigureH);