function varargout = Gui_Tracking(varargin)


% Author: Antoine Godin
% godin.antoine@sympatico.ca

% GUI_TRACKING M-file for Gui_Tracking.fig
%      GUI_TRACKING, by itself, creates a new GUI_TRACKING or raises the existing
%      singleton*.
%
%      H = GUI_TRACKING returns the handle to a new GUI_TRACKING or the handle to
%      the existing singleton*.
%
%      GUI_TRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRACKING.M with the given input arguments.
%
%      GUI_TRACKING('Property','Value',...) creates a new GUI_TRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gui_Tracking_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_Tracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Gui_Tracking

% Last Modified by GUIDE v2.5 12-Jul-2007 18:21:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_Tracking_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_Tracking_OutputFcn, ...
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
                    
% End initialization code - DO NOT EDIT


% --- Executes just before Gui_Tracking is made visible.
function Gui_Tracking_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gui_Tracking (see VARARGIN)

% Choose default command line output for Gui_Tracking
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

            axes(handles.FigureTracks);         axis off
            axes(handles.FigureSeg)   ;         axis off             
            hereIni = cd;
            cd ..
                here = cd;
            cd(hereIni);
            
            if length(varargin) > 0
                handles.FolderName = varargin{1};
                set(handles.ImagesFolder,'String',[handles.FolderName]);           
            else
                handles.FolderName = uigetdir([here,'\']);
                set(handles.ImagesFolder,'String',[handles.FolderName,'\']);           
            end            
            if isa(handles.FolderName,'char') == 1
                if length(dir([ handles.FolderName,'\Settings.ini'])) == 0
                    h = warndlg('Settings.ini does not exist in given folder. Defaults Settings Copied!','Settings.ini?')
                    copyfile('SettingsDefaults.ini',[handles.FolderName,'\Settings.ini'])
                end     
                STRfolder = dir([handles.FolderName,'\*.tif']);        
                if length(STRfolder) > 0
                    for it = 1:length(STRfolder)
                        name = STRfolder(it).name;
                        ima(it) = str2num(name(end-7:end-4));
                    end 
                    NameImages = name(1:length(name)-8);
                    settingsRead=inifile([handles.FolderName,'\Settings.ini'],'readall');    
                            settingsWrite = settingsRead;
                            settingsWrite{1,4} = NameImages;
                            settingsWrite{2,4} = num2str(max(ima)-min(ima)+1);                            
                            inifile([handles.FolderName,'\Settings.ini'],'write',settingsWrite,'plain')
                            
                            set(handles.maxDisplacement,'String',settingsWrite{13,4});
                            set(handles.minLength,'String',settingsWrite{14,4});
                            set(handles.lostFrames ,'String',settingsWrite{15,4});         
                            set(handles.ImageSlider,'Min',min(ima));
                            set(handles.ImageSlider,'Max',max(ima));
                            set(handles.ImageSlider,'Value',min(ima));
                            set(handles.ImageNumber,'String',num2str(min(ima)));
                            set(handles.ImageSlider,'SliderStep',[1/(max(ima)-min(ima)),1/(max(ima)-min(ima))]);
                            if length(settingsRead{17,4}) > 0
                                StringFiloOUTTEMP = str2num(settingsRead{17,4});
                                [a,b] = sort(StringFiloOUTTEMP(:,1));
                                StringFiloOUT = num2str(StringFiloOUTTEMP(b,:));   
                                set(handles.FilopodiaTXT,'String',StringFiloOUT);
                            end
                else
                    h = warndlg('No Images in given Folder!','No Images in folder?')
                end                                 
            end
            guidata(hObject,handles);
% UIWAIT makes Gui_Tracking wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gui_Tracking_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ImagesFolder_Callback(hObject, eventdata, handles)
% hObject    handle to ImagesFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImagesFolder as text
%        str2double(get(hObject,'String')) returns contents of ImagesFolder as a double
            handles.FolderName = get(handles.ImagesFolder,'String');        
            if isa(handles.FolderName,'char') == 1
                if length(dir([ handles.FolderName,'\Settings.ini'])) == 0
                    h = warndlg('Settings.ini does not exist in given folder. Defaults Settings Copied!','Settings.ini?')
                    copyfile('SettingsDefaults.ini',[handles.FolderName,'\Settings.ini'])
                end     
                STRfolder = dir([handles.FolderName,'\*.tif']);        
                if length(STRfolder) > 0
                    for it = 1:length(STRfolder)
                        name = STRfolder(it).name;
                        ima(it) = str2num(name(end-7:end-4));
                    end 
                    NameImages = name(1:length(name)-8);
                            settingsRead=inifile([handles.FolderName,'\Settings.ini'],'readall');    
                            settingsWrite = settingsRead;
                            settingsWrite{1,4} = NameImages;
                            settingsWrite{2,4} = num2str(max(ima)-min(ima)+1);
                            inifile([handles.FolderName,'\Settings.ini'],'write',settingsWrite,'plain')   
                        set(handles.maxDisplacement,'String',settingsWrite{13,4});
                        set(handles.minLength,'String',settingsWrite{14,4});
                        set(handles.lostFrames ,'String',settingsWrite{15,4});             
                            set(handles.ImageSlider,'Min',min(ima));
                            set(handles.ImageSlider,'Max',max(ima));
                            set(handles.ImageSlider,'Value',min(ima));
                            set(handles.ImageNumber,'String',num2str(min(ima)));     
                else
                    h = warndlg('No Images in given Folder!','No Images in folder?')
                end
           end
            

% --- Executes during object creation, after setting all properties.
function ImagesFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImagesFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lostFrames_Callback(hObject, eventdata, handles)
% hObject    handle to lostFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lostFrames as text
%        str2double(get(hObject,'String')) returns contents of lostFrames as a double


% --- Executes during object creation, after setting all properties.
function lostFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lostFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minLength_Callback(hObject, eventdata, handles)
% hObject    handle to minLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLength as text
%        str2double(get(hObject,'String')) returns contents of minLength as a double


% --- Executes during object creation, after setting all properties.
function minLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxDisplacement_Callback(hObject, eventdata, handles)
% hObject    handle to maxDisplacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxDisplacement as text
%        str2double(get(hObject,'String')) returns contents of maxDisplacement as a double


% --- Executes during object creation, after setting all properties.
function maxDisplacement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxDisplacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImageNumber_Callback(hObject, eventdata, handles)
% hObject    handle to ImageNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageNumber as text
%        str2double(get(hObject,'String')) returns contents of ImageNumber as a double
            imagesFolder = get(handles.ImagesFolder,'String');
            frameToTest = floor(str2num(get(handles.ImageNumber,'String')));
            settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');    
            STRfolder = dir([imagesFolder,'Tracks\*.jpg']);        
            if size(STRfolder,1) > 0
                minIma = 0;                
                maxIma = size(STRfolder,1)-1;
                set(handles.ImageSlider,'Min',minIma);
                set(handles.ImageSlider,'Max',maxIma);
                set(handles.ImageSlider,'SliderStep',[1/(maxIma-minIma),1/(maxIma-minIma)]);
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
                ThisImageTrack   = imread([imagesFolder,'Tracks\'  ,settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);
                ThisImageSegment = imread([imagesFolder,'Segments\',settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);

                axes(handles.FigureTracks);         image(ThisImageTrack); axis off;
                axes(handles.FigureSeg)   ;         image(ThisImageSegment); axis off;           
            else
                h = warndlg('No Tracks and Segments!','Tracks and Segments?')
            end

% --- Executes during object creation, after setting all properties.
function ImageNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            imagesFolder = get(handles.ImagesFolder,'String');
            settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');    
            STRfolder = dir([imagesFolder,'Tracks\*.jpg']);        
            if size(STRfolder,1) > 0
                minIma = 0;                
                maxIma = size(STRfolder,1)-1;
                set(handles.ImageSlider,'Min',minIma);
                set(handles.ImageSlider,'Max',maxIma);
                set(handles.ImageSlider,'SliderStep',[1/(maxIma-minIma),1/(maxIma-minIma)]);
                frameToTest = get(handles.ImageSlider,'Value');
                set(handles.ImageSlider,'Value',floor(frameToTest));   
                set(handles.ImageNumber,'String',floor(frameToTest));   
            
                ThisImageTrack   = imread([imagesFolder,'Tracks\'  ,settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);
                ThisImageSegment = imread([imagesFolder,'Segments\',settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);
                axes(handles.FigureTracks);         image(ThisImageTrack); axis off;
                axes(handles.FigureSeg)   ;         image(ThisImageSegment); axis off;           
            else
                h = warndlg('No Tracks and Segments!','Tracks and Segments?')
            end

% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in MakeTrackSegments.
function MakeTrackSegments_Callback(hObject, eventdata, handles)
% hObject    handle to MakeTrackSegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    imagesFolder = get(handles.ImagesFolder,'String');
    settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');    
            settingsWrite=settingsRead;
            settingsWrite{13,4} = get(handles.maxDisplacement,'String');
            settingsWrite{14,4} = get(handles.minLength,'String');
            settingsWrite{15,4} = get(handles.lostFrames,'String');
            inifile([imagesFolder 'Settings.ini'],'write',settingsWrite,'plain')
            
    trackSkeletons
    makeSegmentsMovie
    ImageNumber = 1;
    ThisImageTrack   = imread([imagesFolder,'Tracks\',settingsRead{1,4},num2str(ImageNumber-1, '%04.f'),'.jpg']);
    ThisImageSegment = imread([imagesFolder,'Segments\',settingsRead{1,4},num2str(ImageNumber-1, '%04.f'),'.jpg']);
    
            axes(handles.FigureTracks);         image(ThisImageTrack); axis off;
            axes(handles.FigureSeg)   ;         image(ThisImageSegment); axis off;           
                STRfolder = dir([imagesFolder,'Tracks\*.jpg']);        
                minIma = 0                
                maxIma = size(STRfolder,1)-1;
                set(handles.ImageSlider,'Min',minIma);
                set(handles.ImageSlider,'Max',maxIma);
                set(handles.ImageSlider,'SliderStep',[1/(maxIma-minIma),1/(maxIma-minIma)]);
                set(handles.ImageSlider,'Value',minIma);
                set(handles.ImageNumber,'String',num2str(minIma));
                    
    handles.ImageSize = size(ThisImageTrack);
    guidata(hObject,handles)

% --- Executes on button press in ChoseFilopodia.
function ChoseFilopodia_Callback(hObject, eventdata, handles)
% hObject    handle to ChoseFilopodia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
           
    imagesFolder = get(handles.ImagesFolder,'String');
    FullInfo = load([imagesFolder,'Tracks\Fullinfo.dat']);
    frameToTest = floor(str2num(get(handles.ImageNumber,'String')));
    
    PoFullInfo = find(FullInfo(:,9) == frameToTest+1);
    StFullInfo = FullInfo(PoFullInfo,10)
    [PoFilopodia,ok] = listdlg('ListString',num2str(StFullInfo),'SelectionMode','Single');
    if ok == 1
        Filopodia =   StFullInfo(PoFilopodia);
        set(handles.ImageSlider,'Value',floor(frameToTest));   
        set(handles.ImageNumber,'String',num2str(floor(frameToTest)));   
        settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');    
        ThisImageTrack   = imread([imagesFolder,'Tracks\'  ,settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);
        ThisImageSegment = imread([imagesFolder,'Segments\',settingsRead{1,4},num2str(frameToTest, '%04.f'),'.jpg']);
                axes(handles.FigureTracks);         image(ThisImageTrack); axis off;
                axes(handles.FigureSeg)   ;         image(ThisImageSegment); axis off;     
                handles.ImageSize = size(ThisImageTrack);
        axes(handles.FigureSeg); 
        ValuesGood = 0;
        while ValuesGood == 0
            axes(handles.FigureSeg);
            [X,Y] = ginput(1)
            if (round(X) > 0 & round(X) < handles.ImageSize(2)) & (round(Y) > 0 & round(Y) < handles.ImageSize(1))
                ValuesGood = 1;
            end
        end
        [Filopodia,round(X),round(Y)]
        StringFiloOUT = str2num(get(handles.FilopodiaTXT,'String'))
        if size(StringFiloOUT,1)*size(StringFiloOUT,2) > 0
            switch length(find(StringFiloOUT(:,1)==Filopodia))
                case 0
                    StringFiloOUT(size(StringFiloOUT,1)+1,:) = [Filopodia,round(X),round(Y)];
                case 1
                    StringFiloOUT(find(StringFiloOUT(1:size(StringFiloOUT,1),1)==Filopodia),:) = [Filopodia,round(X),round(Y)];
            end        
        else
            StringFiloOUT =  [Filopodia,round(X),round(Y)];
        end

        [a,b] = sort(StringFiloOUT(:,1));
        StringFiloOUT = num2str(StringFiloOUT(b,:));   
        set(handles.FilopodiaTXT,'String',StringFiloOUT);
        settingsRead = inifile([imagesFolder 'Settings.ini'],'readall');    
        settingsWrite=settingsRead;
        settingsWrite{17,4} = str2num(StringFiloOUT);
        inifile([imagesFolder 'Settings.ini'],'write',settingsWrite,'plain')
    end    
    
% --- Executes on button press in ClearFilopodiaTXT.
function ClearFilopodiaTXT_Callback(hObject, eventdata, handles)
% hObject    handle to ClearFilopodiaTXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    imagesFolder = get(handles.ImagesFolder,'String');
    set(handles.FilopodiaTXT,'String','');
    settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');    
    settingsRead{17,4} = '';
    inifile([imagesFolder 'Settings.ini'],'write',settingsRead,'plain')



function FilopodiaTXT_Callback(hObject, eventdata, handles)
% hObject    handle to FilopodiaTXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilopodiaTXT as text
%        str2double(get(hObject,'String')) returns contents of FilopodiaTXT as a double


% --- Executes during object creation, after setting all properties.
function FilopodiaTXT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilopodiaTXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MakeSegmentsWithArbitraryPoints.
function MakeSegmentsWithArbitraryPoints_Callback(hObject, eventdata, handles)
% hObject    handle to MakeSegmentsWithArbitraryPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    imagesFolder = get(handles.ImagesFolder,'String');
    makeSegmentsWithArbitraryPoints
% --- Executes on button press in FollowSkeleton.
function FollowSkeleton_Callback(hObject, eventdata, handles)
% hObject    handle to FollowSkeleton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    imagesFolder = get(handles.ImagesFolder,'String');
    followSkeletonAnalysis_Point_temp
    


% --- Executes on button press in SeeResults.
function SeeResults_Callback(hObject, eventdata, handles)
% hObject    handle to SeeResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
    Gui_ShowResults(get(handles.ImagesFolder,'String'))
