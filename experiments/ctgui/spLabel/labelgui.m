function varargout = labelgui(varargin)
% LABELGUI Application M-file for labelgui.fig
%    FIG = LABELGUI launch labelgui GUI.
%    LABELGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 21-Oct-2003 14:27:29

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
    
    %initialize colormap
    cMap=hsv(64); %minimal length of cMap: 2*maxColor
    try
        SetUserData(fig,cMap,0);
    catch
        disp('warning: there is already a labelgui open!')
    end
    
    %---------IMPORTANT:CHANGE IF NEW VERSION OF IDLIST NECESSARY TO RUN LABELGUI-------------%
    latestChange='12-Jun-2003';
    SetUserData(fig,latestChange,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
    end
    
end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------Max Intensity radio button------------------------------------------
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton1.
set(handles.radiobutton2,'Value',0);
set(handles.slider4,'Enable','Off'); % stack slider 
set(handles.text2,'Enable','Off');    % stack postion  text
set(handles.text5,'Enable','Off');    % stack sections text
set(h,'Value',1);

imgFigureH=GetUserData(handles.labelgui,'currentWindow');
if isempty(imgFigureH)
    return;
end;

mv=GetUserData(imgFigureH,'mv');

% set slice mode
timesliderH=handles.slider3;
tp=get(timesliderH,'Value');

stacksliderH=handles.slider4;
sp=get(stacksliderH,'Value');

%display z-slice
mIntp=max(mv(:,:,:,1,tp),[],3);
figure(imgFigureH);
cla;

%check wheter we use true brightness or whether the max intensity of
%this image is 1
useTrueBrightness = strcmp(get(handles.label_trueBrightness,'Checked'),'on');
if ~useTrueBrightness
    mIntp = mIntp/max(mIntp(:));
end

imshow(mIntp);
plotlabels(tp);

%if xz/yz-viewer is open: refresh
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    label_refreshXZYZ;
end

%if movieDataFigure is open: refresh
dataFigH = GetUserData(imgFigureH,'CurrentMovieData');    
if ishandle(dataFigH)
    label_showMovieData([],[],handles,1);
end

% ---------------------------Slice mode radio button-----------------------------------------
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton2.
set(handles.radiobutton1,'Value',0);
set(handles.slider4,'Enable','On'); % stack slider 
set(handles.text2,'Enable','On');    % stack postion  text
set(handles.text5,'Enable','On');    % stack sections text
set(h,'Value',1);

imgFigureH=GetUserData(handles.labelgui,'currentWindow');
if isempty(imgFigureH)
    return;
end;

mv=GetUserData(imgFigureH,'mv');

% set slice mode
timesliderH=handles.slider3;
tp=get(timesliderH,'Value');

stackliderH=handles.slider4;
val=get(stackliderH,'Value');
%display z-slice
imgH=findall(imgFigureH,'Type','image');
figure(imgFigureH);
cla;

img = mv(:,:,val,1,tp);
%check wheter we use true brightness or whether the max intensity of
%this image is 1
useTrueBrightness = strcmp(get(handles.label_trueBrightness,'Checked'),'on');
if ~useTrueBrightness
    img = img/max(img(:));
end

imshow(img);
plotlabels(tp,val);

%if xz/yz-viewer is open: refresh
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    label_refreshXZYZ;
end

%if movieDataFigure is open: refresh
dataFigH = GetUserData(imgFigureH,'CurrentMovieData');     
if ishandle(dataFigH)
    label_showMovieData([],[],handles,1);
end

% ------------------------timepoint slider--------------------------------------------
function varargout = slider3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider3.

imgFigureH=GetUserData(handles.labelgui,'currentWindow');
if isempty(imgFigureH)
    return;
end;

mv=GetUserData(imgFigureH,'mv');

val=round(get(h,'Value'));
set(h,'Value',val);
timetextH=handles.text1;
set(timetextH,'String',num2str(val));

%check mode
if get(handles.radiobutton1,'Value')
    mIntp=max(mv(:,:,:,1,val),[],3);
    
    %check wheter we use true brightness or whether the max intensity of
    %this image is 1
    useTrueBrightness = strcmp(get(handles.label_trueBrightness,'Checked'),'on');
    if ~useTrueBrightness
        mIntp = mIntp/max(mIntp(:));
    end
    
    figure(imgFigureH);
    cla;
    %uiViewPanelShowImg(img(:,:,sp,1,val), 0, imgFigureH);
    imshow(mIntp);
    plotlabels(val);
    
else
    
    stacksliderH=handles.slider4;
    sp=get(stacksliderH,'Value');
    
    %display z-slice
    %imgH=findall(imgFigureH,'Type','image');
    %set(imgH,'CData',img(:,:,sp,1,val));
    
    figure(imgFigureH);
    cla;
    
    img = mv(:,:,sp,1,val);
    
    %check wheter we use true brightness or whether the max intensity of
    %this image is 1
    useTrueBrightness = strcmp(get(handles.label_trueBrightness,'Checked'),'on');
    if ~useTrueBrightness
        img = img/max(img(:));
    end
    
    imshow(img);
    plotlabels(val,sp);
    
end;

%if connectGUI is open: reopen it with a new time
cguiH=findall(0,'Tag','connectGUI');
if ishandle(cguiH)
    close(cguiH);
    connectGUI;
end

%if view3D is open & labelgui sets the time, adjust it
view3DH=GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    currTimeByLabelgui=get(findall(0,'Tag','view3D_setTime_PD'),'Value');
    if currTimeByLabelgui==2
        view3D_refresh;
    end
end

%if xz/yz-viewer is open: refresh
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    label_refreshXZYZ;
end

%if movieDataFigure is open: refresh
dataFigH = GetUserData(imgFigureH,'CurrentMovieData');   
if ishandle(dataFigH)
    label_showMovieData([],[],handles,1);
end

% ----------------------------section slider----------------------------------------
function varargout = slider4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider4.

imgFigureH=GetUserData(handles.labelgui,'currentWindow');
if isempty(imgFigureH)
    return;
end;
mv=GetUserData(imgFigureH,'mv');

val=round(get(gcbo,'Value'));
set(gcbo,'Value',val);
stacktextH=handles.text2;
timesliderH=handles.slider3;
tp=get(timesliderH,'Value');

set(stacktextH,'String',num2str(val));

%display z-slice
imgH=findall(imgFigureH,'Type','image');
figure(imgFigureH);

img = mv(:,:,val,1,tp);
    
    %check wheter we use true brightness or whether the max intensity of
    %this image is 1
    useTrueBrightness = strcmp(get(handles.label_trueBrightness,'Checked'),'on');
    if ~useTrueBrightness
        img = img/max(img(:));
    end

set(imgH,'CData',img);

plotlabels(tp,val);

%if xz/yz-viewer is open: refresh
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    label_refreshXZYZ;
end

%if movieDataFigure is open: refresh
dataFigH = GetUserData(imgFigureH,'CurrentMovieData');   
if ishandle(dataFigH)
    label_showMovieData([],[],handles,1);
end

% ----------------------------plot labels----------------------------------------
function plotlabels(timepoint,section);

%read userData
labelGuiH=findall(0,'Tag','labelgui');
imgFigureH=GetUserData(labelGuiH,'currentWindow');
idlist=GetUserData(imgFigureH,'idlist');
if isempty(idlist)
    return %no idlist loaded -> no plotlabels
end
cMap=GetUserData(labelGuiH,'cMap');
cMapFact=size(cMap,1)/idlist(1).stats.maxColor;

%get constants
dataProperties=GetUserData(imgFigureH,'dataProperties');
if isempty(dataProperties)
    error('no dataProperties found');
    %     load_sp_constants;
    %     global PIXELSIZE_XY PIXELSIZE_Z;
else
    PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
    PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;
end

%for whatever reason the gca are no longer found running the soft from
%runctbatch. hence make everything visible for a moment
set(imgFigureH,'HandleVisibility','on');

if ~isempty(idlist(timepoint).linklist)
    for i=1:max(idlist(timepoint).linklist(:,2))
        
        %find colors belonging to same spot
        rowIdx=find(idlist(timepoint).linklist(:,2)==i);
        
        %read and change coordinates
        coord=idlist(timepoint).linklist(rowIdx(1),9:11)./[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
        
        %list of colors in one spot
        cList=idlist(timepoint).linklist(rowIdx,4);
        
        if (nargin==1 | round(coord(3))==section)
            if length(cList)==1
                % plot position
                pH=plot(coord(1),coord(2),'o','MarkerSize',4);
                labelColor(1,:)=cMap(cList*cMapFact,:);
                set(pH,'MarkerEdgeColor',labelColor);
                type=char(idlist(1).stats.labelcolor(log2(cList)+1));
            else
                %plot all spots making up this one
                for j=length(cList):-1:1
                    pH=plot(coord(1),coord(2),'o','MarkerSize',4+4*(j-1));
                    labelColor(j,:)=cMap(cList(j)*cMapFact,:);
                    set(pH,'MarkerEdgeColor',labelColor(j,:));
                end
                labelColor=mean(labelColor,1);
                type='fus';
            end
            
            % write label
            txtH=text(coord(1)+5,coord(2),type,'Color',labelColor);
            %set label2spot number
            set(txtH,'UserData',i);
            % popup CB
            set(txtH,'ButtonDownFcn','popuphere');    
        end;
    end;  
end;

%and hide the figure again
set(imgFigureH,'HandleVisibility','callback');

% make gui top most
figure(labelGuiH);

%------------------------------------------------------------------------------------

function refresh

imgFigureH=GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    return;
end;
figure(imgFigureH);

% remove old labels
lineH=findall(imgFigureH,'Type','line');
txtH=findall(imgFigureH,'Type','text');
delete([txtH;lineH]);

cfigH=findall(0,'Tag','labelgui');
handles=guidata(cfigH);
tp=get(handles.slider3,'Value');
if(get(handles.radiobutton2,'Value'))
    slice=get(handles.slider4,'Value');
    plotlabels(tp,slice);
else
    plotlabels(tp);
end;

%make the right guis visible
view3DH=GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    set(view3DH,'Visible','on')
    view3DfigH = findall(0,'Tag','view3Dfig');
    if ishandle(view3DfigH)
        set(view3DfigH,'Visible','off');
    end
    currTimeByLabelgui=get(findall(0,'Tag','view3D_setTime_PD'),'Value');
    if currTimeByLabelgui==2
        view3D_refresh;
    end
end

%if xz/yz-viewer is open: refresh
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    set(zFigH,'Visible','on');
    label_refreshXZYZ;
end

%if movieDataFigure is open: refresh
dataFigH = GetUserData(imgFigureH,'CurrentMovieData');    
if ishandle(dataFigH)
    label_showMovieData([],[],handles,1);
    set(dataFigH,'Visible','on');
end

%intFigure
figH = GetUserData(imgFigureH,'intFigH');
if ishandle(figH)
    set(figH,'Visible','on');
end

%----------------------------------------------------------------------
function label_chooseWindow_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_chooseWindow_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in label_chooseWindow_PD.
function label_chooseWindow_PD_Callback(hObject, eventdata, handles)

%switch: if load movie, go to label_load
%if load idlist, goto load idlist
pdPos = get(hObject,'Value');
lengthPDList = size(get(hObject, 'String'),1);

switch lengthPDList-pdPos
    %if pdPos==lengthPDList: load movie
    %if pdPos==lengthPDList-1: load idlist
    %else: do nothing
    case 0 %load movie
        label_loadCB;
    case 1
        label_loadslistCB;
    otherwise
        
        %make sure everything about the newly active figure is ok
        %and store the right handle
        handleList = GetUserData(handles.labelgui,'labelPanelHandles');
        
        %check that the user did not just not change anything
        if pdPos == 1 & isempty(handleList)
            return
        end
        
        imgFigureH = handleList(pdPos);
        SetUserData(handles.labelgui,imgFigureH,1,'currentWindow');
        
        %make all additional figures invisible
        %3d-viewer
        view3DH = findall(0,'Tag','view3Dgui');
        if ishandle(view3DH)
            set(view3DH,'Visible','off');
        end
        view3DfigH = findall(0,'Tag','view3Dfig');
        if ishandle(view3DfigH)
            set(view3DfigH,'Visible','off');
        end
        
%         %xz/yz-viewer
%         zFigH = findall(0,'Tag','XZYZFigure');
%         if ishandle(zFigH)
%             set(zFigH,'Visible','off');
%         end
        
        %movieDataFigure
        dataFigH = findall(0,'Name','CurrentMovieData');    
        if ishandle(dataFigH)
            set(dataFigH,'Visible','off');
        end
        
        %intFigure
        figH = findall(0,'Tag','intFig');
        if ishandle(figH)
            set(figH,'Visible','off');
        end
        
        %refresh all: hence call tpslider/stackslider/refresh
        %timepoint slider
        labelgui('slider3_Callback',handles.slider3,[],handles);
        %stack slider
        if strcmp(get(handles.slider4,'Enable'),'On')
            labelgui('slider4_Callback',handles.slider4,[],handles);
        end
        
        labelgui('refresh');
        
end


% --------------------------------------------------------------------
function label_trueBrightness_Callback(hObject, eventdata, handles)
% hObject    handle to label_trueBrightness_CB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%unlike a checkbox, a checkMenu has to be (un)checked manually
isChecked = strcmp(get(hObject,'Checked'),'on');
if isChecked
    set(hObject,'Checked','off')
else
    set(hObject,'Checked','on')
end

%whether this is being checked or not - the figure has to be updated (if there is one)
imgFigureH=GetUserData(handles.labelgui,'currentWindow');
if isempty(imgFigureH)
    return;
end;

%do not update with refresh, but with the timepointSlider, because the
%image has to be calculated again
labelgui('slider3_Callback',handles.slider3, eventdata, handles);


