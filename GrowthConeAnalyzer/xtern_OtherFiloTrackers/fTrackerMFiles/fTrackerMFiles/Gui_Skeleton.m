function varargout = Gui_Skeleton(varargin)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

    % GUI_SKELETON M-file for Gui_Skeleton.fig
    %      GUI_SKELETON, by itself, creates a new GUI_SKELETON or raises the existing
    %      singleton*.
    %
    %      H = GUI_SKELETON returns the handle to a new GUI_SKELETON or the handle to
    %      the existing singleton*.
    %
    %      GUI_SKELETON('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in GUI_SKELETON.M with the given input arguments.
    %
    %      GUI_SKELETON('Property','Value',...) creates a new GUI_SKELETON or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before Gui_Skeleton_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to Gui_Skeleton_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Copyright 2002-2003 The MathWorks, Inc.

    % Edit the above text to modify the response to help Gui_Skeletton

    % Last Modified by GUIDE v2.5 09-Jul-2007 13:23:16

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Gui_Skeleton_OpeningFcn, ...
                       'gui_OutputFcn',  @Gui_Skeleton_OutputFcn, ...
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

    % --- Executes just before Gui_Skeletton is made visible.
function Gui_Skeleton_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to Gui_Skeletton (see VARARGIN)

    % Choose default command line output for Gui_Skeletton
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    % UIWAIT makes Gui_Skeleton wait for user response (see UIRESUME)
    % uiwait(handles.Gui_Skeleton);

        axes(handles.Figure1)    ;             axis off  
        axes(handles.Figure2)    ;             axis off  
        axes(handles.Figure3)    ;             axis off  
        axes(handles.Figure4)    ;             axis off  
        axes(handles.Figure5)    ;             axis off  
        axes(handles.Figure6)    ;             axis off   
        axes(handles.Figure7)    ;             axis off  
        axes(handles.Figure8)    ;             axis off  
        hereIni = cd;
        cd ..
        here = cd;
        cd(hereIni);
        handles.FolderName = uigetdir([here,'\']);
                
        set(handles.imagesFolder,'String',[handles.FolderName,'\']);
        imagesFolder = handles.FolderName;
        if isa(handles.FolderName,'char') == 1
            if length(dir([ handles.FolderName,'\Settings.ini'])) == 0
                h = warndlg('Settings.ini does not exist in given folder. Defaults Settings Copied!','Settings.ini?')
                copyfile('SettingsDefaults.ini',[handles.FolderName,'\Settings.ini'])
            end
        end
                
        STRfolder = dir([handles.FolderName,'\*.tif']);        
        if length(STRfolder) > 0
            try
                for it = 1:length(STRfolder)
                    name = STRfolder(it).name;
                    ima(it) = str2num(name(end-7:end-4));
                end 
                NameImages = name(1:length(name)-8);
            catch
                for it = 1:length(STRfolder)
                    name = STRfolder(it).name;
                    ima(it) = str2num(name(end-6:end-4));
                end 
                NameImages = name(1:length(name)-7);
                for it = min(ima):max(ima)
                    if min(ima) == 0                         
                        [status,message,messageid] = copyfile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'], ...
                                                              [handles.FolderName '\' NameImages num2str(it, '%04.f') '.tif'],'f');
                    else
                        [status,message,messageid] = copyfile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'], ...
                                                              [handles.FolderName '\' NameImages num2str(it-min(ima), '%04.f') '.tif'],'f');
                    end
                        [i j k] = mkdir([handles.FolderName '\OldFiles']);
                        movefile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'],[handles.FolderName '\OldFiles\' NameImages num2str(it, '%03.f') '.tif'])
%                         delete([Foldertotest NameImages num2str(it, '%03.f') '.tif']);       
                end
                ima = ima-min(ima);
            end
            set(handles.FrameToTestSlider,'Min',min(ima));
            set(handles.FrameToTestSlider,'Max',max(ima));
            set(handles.FrameToTestSlider,'SliderStep',[1/(max(ima)-min(ima)),1/(max(ima)-min(ima))]);

            settingsRead=inifile([handles.FolderName,'\Settings.ini'],'readall');    
                    settingsWrite = settingsRead;
                    settingsWrite{1,4} = NameImages;
                    settingsWrite{2,4} = num2str(max(ima)-min(ima)+1);
                    inifile([handles.FolderName,'\Settings.ini'],'write',settingsWrite,'plain')    
               
                set(handles.ThresholdFactor,'String',settingsWrite{3,4});
                set(handles.PixelsDilate   ,'String',settingsWrite{4,4});
                set(handles.DiskSize       ,'String',settingsWrite{6,4});
                set(handles.CloseSize      ,'String',settingsWrite{7,4});
                set(handles.EdgeStrongness ,'String',settingsWrite{9,4});
                
                set(handles.SkelFill       ,'Value',str2num(settingsWrite{5,4}));
                set(handles.IntensityBased ,'Value',str2num(settingsWrite{8,4}));
                if get(handles.IntensityBased,'Value') == 1
                    set(handles.ThresholdFactor,'Enable','on')
                    set(handles.DiskSize,'Enable','on')
                else
                    set(handles.ThresholdFactor,'Enable','off')
                    set(handles.DiskSize,'Enable','off')
                end                
        else
            h = warndlg('No Images in given Folder!','No Images in folder?')
        end
        guidata(hObject,handles);

    % --- Outputs from this function are returned to the command line.
function varargout = Gui_Skeleton_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

function imagesFolder_Callback(hObject, eventdata, handles)
    % hObject    handle to imagesFolder (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of imagesFolder as text
    %        str2double(get(hObject,'String')) returns contents of imagesFolder as a double
    imagesFolder_CreateFcn(hObject, eventdata, handles)
    
    % --- Executes during object creation, after setting all properties.
function imagesFolder_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to imagesFolder (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
        if isa(handles,'struct')==1
            Err = Check_Inputs(hObject, eventdata, handles);        
            if Err == 0        
                imagesFolder = get(handles.imagesFolder,'String');
                handles.FolderName = get(handles.imagesFolder,'String');        
                if isa(handles.FolderName,'char') == 1
                    if length(dir([ handles.FolderName,'\Settings.ini'])) == 0
                        h = warndlg('Settings.ini does not exist in given folder. Defaults Settings Copied!','Settings.ini?')
                        copyfile('SettingsDefaults.ini',[handles.FolderName,'\Settings.ini'])
                    end
                end
       
                STRfolder = dir([handles.FolderName,'\*.tif']);        
                try
                    for it = 1:length(STRfolder)
                        name = STRfolder(it).name;
                        ima(it) = str2num(name(end-7:end-4));
                    end 
                    NameImages = name(1:length(name)-8);
                catch
                    for it = 1:length(STRfolder)
                        name = STRfolder(it).name;
                        ima(it) = str2num(name(end-6:end-4));
                    end 
                    NameImages = name(1:length(name)-7);
                    for it = min(ima):max(ima)
                        if min(ima) == 0                         
                            [status,message,messageid] = copyfile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'], ...
                                                                  [handles.FolderName '\' NameImages num2str(it, '%04.f') '.tif'],'f');
                        else
                            [status,message,messageid] = copyfile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'], ...
                                                                  [handles.FolderName '\' NameImages num2str(it-min(ima), '%04.f') '.tif'],'f');
                        end
                            [i j k] = mkdir([handles.FolderName '\OldFiles']);
                            movefile([handles.FolderName '\' NameImages num2str(it, '%03.f') '.tif'],[handles.FolderName '\OldFiles\' NameImages num2str(it, '%03.f') '.tif'])
                    end
                    ima = ima-min(ima);
                end
                    set(handles.FrameToTestSlider,'Min',min(ima));
                    set(handles.FrameToTestSlider,'Max',max(ima));
                    set(handles.FrameToTestSlider,'SliderStep',[1/(max(ima)-min(ima)),1/(max(ima)-min(ima))]);

                    settingsRead=inifile([handles.FolderName,'\Settings.ini'],'readall');    
                            settingsWrite = settingsRead;
                            settingsWrite{1,4} = NameImages;
                            settingsWrite{2,4} = num2str(max(ima)-min(ima)+1);
                            inifile([handles.FolderName,'\Settings.ini'],'write',settingsWrite,'plain')    
                            
                    set(handles.ThresholdFactor,'String',settingsWrite{3,4});
                    set(handles.PixelsDilate,'String',settingsWrite{4,4});
                    set(handles.DiskSize ,'String',settingsWrite{6,4});
                    set(handles.CloseSize,'String',settingsWrite{7,4});
                    
                    set(handles.SkelFill,'Value',str2num(settingsWrite{5,4}));
                    set(handles.IntensityBased,'Value',str2num(settingsWrite{8,4}));       
                    if get(handles.IntensityBased,'Value') == 1
                        set(handles.ThresholdFactor,'Enable','on')
                    else
                        set(handles.ThresholdFactor,'Enable','off')
                    end
                else
                    h = warndlg('No Images in given Folder!','No Images in folder?')
                end
            end

    % --- Executes on button press in AnalyzeFrame.
    function AnalyzeFrame_Callback(hObject, eventdata, handles)
    % hObject    handle to AnalyzeFrame (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

        RunSingleImage_CreateFnc(hObject, eventdata, handles)


                function ThresholdFactor_Callback(hObject, eventdata, handles)
                % hObject    handle to ThresholdFactor (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of ThresholdFactor as text
                %        str2double(get(hObject,'String')) returns contents of ThresholdFactor as a double


                % --- Executes during object creation, after setting all properties.
                function ThresholdFactor_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to ThresholdFactor (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end



                function PixelsDilate_Callback(hObject, eventdata, handles)
                % hObject    handle to PixelsDilate (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of PixelsDilate as text
                %        str2double(get(hObject,'String')) returns contents of PixelsDilate as a double


                % --- Executes during object creation, after setting all properties.
                function PixelsDilate_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to PixelsDilate (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end



                function SkelFill_Callback(hObject, eventdata, handles)
                % hObject    handle to SkelFill (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of SkelFill as text
                %        str2double(get(hObject,'String')) returns contents of SkelFill as a double


                % --- Executes during object creation, after setting all properties.
                function SkelFill_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to SkelFill (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end



                function DiskSize_Callback(hObject, eventdata, handles)
                % hObject    handle to DiskSize (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of DiskSize as text
                %        str2double(get(hObject,'String')) returns contents of DiskSize as a double


                % --- Executes during object creation, after setting all properties.
                function DiskSize_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to DiskSize (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end



                function CloseSize_Callback(hObject, eventdata, handles)
                % hObject    handle to CloseSize (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of CloseSize as text
                %        str2double(get(hObject,'String')) returns contents of CloseSize as a double


                % --- Executes during object creation, after setting all properties.
                function CloseSize_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to CloseSize (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end



                function IntensityBased_Callback(hObject, eventdata, handles)
                % hObject    handle to IntensityBased (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Hints: get(hObject,'String') returns contents of IntensityBased as text
                %        str2double(get(hObject,'String')) returns contents of IntensityBased as a double
                if get(handles.IntensityBased,'Value') == 1
                    set(handles.ThresholdFactor,'Enable','on')     
                    set(handles.DiskSize,'Enable','on')
                    set(handles.EdgeStrongness,'Enable','off')
                else
                    set(handles.ThresholdFactor,'Enable','off')
                    set(handles.DiskSize,'Enable','off')
                    set(handles.EdgeStrongness,'Enable','on')
                end
                    
                % --- Executes during object creation, after setting all properties.
                function IntensityBased_CreateFcn(hObject, eventdata, handles)
                % hObject    handle to IntensityBased (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    empty - handles not created until after all CreateFcns called

                % Hint: edit controls usually have a white background on Windows.
                %       See ISPC and COMPUTER.
                if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                    set(hObject,'BackgroundColor','white');
                end

    function frameToTest_Callback(hObject, eventdata, handles)
        % hObject    handle to frameToTest (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % Hints: get(hObject,'String') returns contents of frameToTest as text
        %        str2double(get(hObject,'String')) returns contents of frameToTest as a double
        Err = Check_Inputs(hObject, eventdata, handles);        
        if Err == 0                  
            frameToTest = str2num(get(handles.frameToTest,'String'));
            set(handles.FrameToTestSlider,'Value',frameToTest);
            imagesFolder = get(handles.imagesFolder,'String');
            settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');           
            fileName = getFromINI(settingsRead, 'fileName');       
            thisFrame=imread([imagesFolder fileName num2str(frameToTest, '%04.f') '.tif']);
                        
            axes(handles.Figure2) 
                    hold off;  image(thisFrame); axis off  
                    set(handles.Figure2TXT,'String','Original Image')
                    set(handles.Figure2TXT,'FontSize',10)           
        end
        % --- Executes during object creation, after setting all properties.
    function frameToTest_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to frameToTest (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: edit controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
    % --- Executes on slider movement.
    function FrameToTestSlider_Callback(hObject, eventdata, handles)
        Err = Check_Inputs(hObject, eventdata, handles);        
        if Err == 0                  
            frameToTest = get(handles.FrameToTestSlider,'Value');
            imagesFolder = get(handles.imagesFolder,'String');

            set(handles.FrameToTestSlider,'Value',floor(frameToTest));   
            set(handles.frameToTest,'String',floor(frameToTest));   
            
            settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');           
            fileName = getFromINI(settingsRead, 'fileName');       
            thisFrame=imread([imagesFolder fileName num2str(frameToTest, '%04.f') '.tif']);
                        
            axes(handles.Figure2) 
                    hold off;  image(thisFrame); axis off  
                    set(handles.Figure2TXT,'String','Original Image')
                    set(handles.Figure2TXT,'FontSize',10)            
        end

    % --- Executes during object creation, after setting all properties.
    function FrameToTestSlider_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to FrameToTestSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB

    % handles    empty - handles not created until after all CreateFcns called

    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


    function RunSingleImage_CreateFnc(hObject, eventdata, handles)
        
        Err = Check_Inputs(hObject, eventdata, handles);
        if Err == 0        
            frameToTest  = str2num(get(handles.frameToTest,'String'));
            imagesFolder = get(handles.imagesFolder,'String');
            
            settingsRead = inifile([imagesFolder 'Settings.ini'],'readall');    
                settingsWrite = settingsRead;
                settingsWrite{3,4} = get(handles.ThresholdFactor,'String');
                settingsWrite{4,4} = get(handles.PixelsDilate,'String');
                settingsWrite{5,4} = get(handles.SkelFill,'Value');
                settingsWrite{6,4} = num2str(round(str2num(get(handles.DiskSize ,'String'))));
                settingsWrite{7,4} = num2str(round(str2num(get(handles.CloseSize,'String'))));
                settingsWrite{8,4} = get(handles.IntensityBased,'Value');              
                settingsWrite{9,4} = num2str(get(handles.EdgeStrongness,'String'));              
                inifile([imagesFolder 'Settings.ini'],'write',settingsWrite,'plain')
            
            singleImageSkeletonAnalysis
            if get(handles.IntensityBased,'Value') == 0
                axes(handles.Figure1)                 
                    AA = zeros(size(imageBackground,1),size(imageBackground,2),3);
                    AA(:,:,1) = double(imageBackground)/double(max(max(imageBackground)));
                    AA(:,:,2) = double(imageBackground)/double(max(max(imageBackground)));
                    AA(:,:,3) = double(imageBackground)/double(max(max(imageBackground)));
                    hold off;  image(AA); axis off  
                    set(handles.Figure1TXT,'String','Image Background')
                    set(handles.Figure1TXT,'FontSize',10)
                axes(handles.Figure2) 
                    hold off;  image(thisFrame); axis off  
                    set(handles.Figure2TXT,'String','Original Image')
                    set(handles.Figure2TXT,'FontSize',10)
                axes(handles.Figure3)     
                    AA = zeros(size(edgesImage,1),size(edgesImage,2),3);
                    AA(:,:,1) = double(edgesImage)/max(max(edgesImage));
                    AA(:,:,2) = double(edgesImage)/max(max(edgesImage));
                    AA(:,:,3) = double(edgesImage)/max(max(edgesImage));
                    hold off;  image(AA); axis off  
                    set(handles.Figure3TXT,'String','Edges Image')               
                    set(handles.Figure3TXT,'FontSize',10)
                axes(handles.Figure4)               
                    AA = zeros(size(closedCone,1),size(closedCone,2),3);
                    AA(:,:,1) = double(closedCone)/max(max(closedCone));
                    AA(:,:,2) = double(closedCone)/max(max(closedCone));
                    AA(:,:,3) = double(closedCone)/max(max(closedCone));
                    hold off;  image(AA);  axis off
                    set(handles.Figure4TXT,'String','Closed Cone')
                    set(handles.Figure4TXT,'FontSize',10)
                axes(handles.Figure5) 
                    AA = zeros(size(dilatedCone,1),size(dilatedCone,2),3);
                    AA(:,:,1) = double(dilatedCone)/double(max(max(dilatedCone)));
                    AA(:,:,2) = double(dilatedCone)/double(max(max(dilatedCone)));
                    AA(:,:,3) = double(dilatedCone)/double(max(max(dilatedCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure5TXT,'String','Dilated Cone')
                    set(handles.Figure5TXT,'FontSize',10)
                axes(handles.Figure6) 
                    AA = zeros(size(singleCone,1),size(singleCone,2),3);
                    AA(:,:,1) = double(singleCone)/double(max(max(singleCone)));
                    AA(:,:,2) = double(singleCone)/double(max(max(singleCone)));
                    AA(:,:,3) = double(singleCone)/double(max(max(singleCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure6TXT,'String','Single Cone')
                    set(handles.Figure6TXT,'FontSize',10)
                axes(handles.Figure7) 
                    AA = zeros(size(filledCone,1),size(filledCone,2),3);
                    AA(:,:,1) = double(filledCone)/double(max(max(filledCone)));
                    AA(:,:,2) = double(filledCone)/double(max(max(filledCone)));
                    AA(:,:,3) = double(filledCone)/double(max(max(filledCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure7TXT,'String','Filled Cone')
                    set(handles.Figure7TXT,'FontSize',10)
                axes(handles.Figure8) 
                    hold off;  imagesc(rgbFrame); axis off             
                    set(handles.Figure8TXT,'String','Final Image')
                    set(handles.Figure8TXT,'FontSize',10)                
            elseif get(handles.IntensityBased,'Value') == 1  
                axes(handles.Figure1) 
                    AA = zeros(size(imageBackground,1),size(imageBackground,2),3);
                    AA(:,:,1) = double(imageBackground)/double(max(max(imageBackground)));
                    AA(:,:,2) = double(imageBackground)/double(max(max(imageBackground)));
                    AA(:,:,3) = double(imageBackground)/double(max(max(imageBackground)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure1TXT,'String','Image Background')
                    set(handles.Figure1TXT,'FontSize',10)
                axes(handles.Figure2) 
                    hold off;  image(thisFrame); axis off  
                    set(handles.Figure2TXT,'String','Original Image')
                    set(handles.Figure2TXT,'FontSize',10)
                axes(handles.Figure3) 
                    AA = zeros(size(binaryImage,1),size(binaryImage,2),3);
                    AA(:,:,1) = double(binaryImage)/double(max(max(binaryImage)));
                    AA(:,:,2) = double(binaryImage)/double(max(max(binaryImage)));
                    AA(:,:,3) = double(binaryImage)/double(max(max(binaryImage)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure3TXT,'String','Binary Image')
                    set(handles.Figure3TXT,'FontSize',10)
                axes(handles.Figure4) 
                    AA = zeros(size(singleCone,1),size(singleCone,2),3);
                    AA(:,:,1) = double(singleCone)/double(max(max(singleCone)));
                    AA(:,:,2) = double(singleCone)/double(max(max(singleCone)));
                    AA(:,:,3) = double(singleCone)/double(max(max(singleCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure4TXT,'String','Single Cone')
                    set(handles.Figure4TXT,'FontSize',10)              
                axes(handles.Figure5) 
                    AA = zeros(size(dilatedCone,1),size(dilatedCone,2),3);
                    AA(:,:,1) = double(dilatedCone)/double(max(max(dilatedCone)));
                    AA(:,:,2) = double(dilatedCone)/double(max(max(dilatedCone)));
                    AA(:,:,3) = double(dilatedCone)/double(max(max(dilatedCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure5TXT,'String','Dilated Cone')
                    set(handles.Figure5TXT,'FontSize',10)
                axes(handles.Figure6) 
                    AA = zeros(size(filteredCone,1),size(filteredCone,2),3);
                    AA(:,:,1) = double(filteredCone)/double(max(max(filteredCone)));
                    AA(:,:,2) = double(filteredCone)/double(max(max(filteredCone)));
                    AA(:,:,3) = double(filteredCone)/double(max(max(filteredCone)));
                    hold off;  image(AA);  axis off
                    set(handles.Figure6TXT,'String','Filtered Cone')
                    set(handles.Figure6TXT,'FontSize',10)
                axes(handles.Figure7) 
                    AA = zeros(size(closedCone,1),size(closedCone,2),3);
                    AA(:,:,1) = double(closedCone)/double(max(max(closedCone)));
                    AA(:,:,2) = double(closedCone)/double(max(max(closedCone)));
                    AA(:,:,3) = double(closedCone)/double(max(max(closedCone)));
                    hold off;  image(AA); axis off
                    set(handles.Figure7TXT,'String','Closed Cone')
                    set(handles.Figure7TXT,'FontSize',10)
                axes(handles.Figure8) 
                    hold off;  imagesc(rgbFrame); axis off             
                    set(handles.Figure8TXT,'String','Final Image')
                    set(handles.Figure8TXT,'FontSize',10)
            end
        end
        'Job Finished'
        guidata(hObject,handles)


% --- Executes on button press in AnalyzeMovie.
function AnalyzeMovie_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    Err = Check_Inputs(hObject, eventdata, handles)
    if Err == 0
        imagesFolder = get(handles.imagesFolder,'String');
        analyzeSkeletonMovie
    end


% --- Executes on button press in OpenGuiTracking.
function OpenGuiTracking_Callback(hObject, eventdata, handles)
% hObject    handle to OpenGuiTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    Gui_Tracking(get(handles.imagesFolder,'String'));


function EdgeStrongness_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function EdgeStrongness_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


