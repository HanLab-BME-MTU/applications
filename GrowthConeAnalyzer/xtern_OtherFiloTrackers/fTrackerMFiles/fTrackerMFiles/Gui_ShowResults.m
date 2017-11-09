function varargout = Gui_ShowResults(varargin)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Gui_ShowResults_OpeningFcn, ...
                       'gui_OutputFcn',  @Gui_ShowResults_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end

function Gui_ShowResults_OpeningFcn(hObject, eventdata, handles, varargin)
    axes(handles.ArbitraryFig);         axis off
    axes(handles.FollowFig)   ;         axis off             
          
    handles.output = hObject;
    handles.FolderName = varargin{1};
    settingsRead=inifile([handles.FolderName 'Settings.ini'],'readall');
    handles.fileName=getFromINI(settingsRead, 'fileName');
    handles.frames=getFromINI(settingsRead, 'frames');
    
                set(handles.ImageSlider,'Min',0);
                set(handles.ImageSlider,'Max',handles.frames-1);
                set(handles.ImageSlider,'SliderStep',[1/(handles.frames-1),1/(handles.frames-1)]);
    frameToTest = 1;
    ThisImageTrack   = imread([handles.FolderName,'Arbitrary\'  ,handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
    ThisImageSegment = imread([handles.FolderName,'FollowSkeleton\',handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
                
    axes(handles.ArbitraryFig);         image(ThisImageTrack); axis off;
    axes(handles.FollowFig)   ;         image(ThisImageSegment); axis off;           
        
    guidata(hObject, handles);
    
function varargout = Gui_ShowResults_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;

function ImageNumber_Callback(hObject, eventdata, handles)
    frameToTest = floor(str2num(get(handles.ImageNumber,'String')));
    minIma = 0;
    maxIma =handles.frames;
    if floor(frameToTest) >= minIma & floor(frameToTest) < maxIma
        set(handles.ImageSlider,'Value',floor(frameToTest));   
        set(handles.ImageNumber,'String',num2str(floor(frameToTest)));   
    elseif floor(frameToTest) < minIma
        frameToTest = minIma;
        set(handles.ImageSlider,'Value',floor(frameToTest));   
        set(handles.ImageNumber,'String',num2str(floor(frameToTest)));   
    else
        frameToTest = maxIma;
        set(handles.ImageSlider,'Value',floor(frameToTest));   
        set(handles.ImageNumber,'String',num2str(floor(frameToTest)));   
    end

    ThisImageTrack   = imread([handles.FolderName,'Arbitrary\'  ,handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
    ThisImageSegment = imread([handles.FolderName,'FollowSkeleton\',handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
                
    axes(handles.ArbitraryFig);         image(ThisImageTrack); axis off;
    axes(handles.FollowFig)   ;         image(ThisImageSegment); axis off;           
    
function ImageNumber_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function ImageSlider_Callback(hObject, eventdata, handles)
    frameToTest = get(handles.ImageSlider,'Value');
    set(handles.ImageSlider,'Value',floor(frameToTest));   
    set(handles.ImageNumber,'String',floor(frameToTest));   

    ThisImageTrack   = imread([handles.FolderName,'Arbitrary\'  ,handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
    ThisImageSegment = imread([handles.FolderName,'FollowSkeleton\',handles.fileName,num2str(frameToTest, '%04.f'),'.jpg']);
                
    axes(handles.ArbitraryFig);         image(ThisImageTrack); axis off;
    axes(handles.FollowFig)   ;         image(ThisImageSegment); axis off;           
    
function ImageSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end