function varargout = fsmGuiMain(varargin)
% fsmGuiMain Application M-file for fsmGuiMain.fig
%    FIG = fsmGuiMain launch fsmGuiMain GUI.
%    fsmGuiMain('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 11-Mar-2004 17:30:09
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  INITIALIZING WINDOW
%  LOADING DEFAULT VALUES (fsmParams), SETTING GUI FIELDS WITH DEFAULT VALUES
%  STORING DEFAULT VALUES INTO defaultFsmParams
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Link fsmParam and defaultFsmParam to the GUI
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load fsmParam structure from fsmParam.mat, if it exists
    pathOfFsmMain=which('fsmMain.m');
    % Get path for fsmParam
    indx=find(pathOfFsmMain==filesep);
    indx=indx(length(indx));
    fsmParamPath=[pathOfFsmMain(1:indx),'fsmParam.mat'];
    if exist(fsmParamPath)==2
        load(fsmParamPath); % Load defaut values stored in the fsm directory
    else
        % Initialize fsmParam with default values
        fsmParam=fsmGetParamDflts;
    end
    
    % Store default values in defaultFsmParam
    defaultFsmParam=fsmParam;
    
    % Link fsmParam to the start button
    set(handles.start,'UserData',fsmParam);
    
    % Link defaultFsmParam to the default button
    set(handles.defaultButton,'UserData',defaultFsmParam);

    % Link z values to the z value edit box
    confidenceProb=[1.15 1.29 1.45 1.645 1.96 2.58];
    set(handles.editZValue,'UserData',confidenceProb);

    % Fill all fields with the values from fsmParam
    fsmGuiWriteParameters(defaultFsmParam,handles);

    % Read parameter experiments from fsmExpParams.txt
    userDir=fsmCenter_getUserSettings;
    if isempty(userDir)
        fsmExpParamPath=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
    else
        fsmExpParamPath=[userDir,filesep,'fsmExpParams.txt'];
        if exist(fsmExpParamPath)~=2
            % No database found in user-defined directory
            % Reverting to default database
            fsmExpParamPath=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
        end
    end
    if exist(fsmExpParamPath)~=2
        uiwait(msgbox('Could not find experiment database! Get fsmExpParams.txt from the repository and restart SpeckTackle.','Error','modal'));
        return
    end
    
    % Fill parameters structure
    fsmExpParam=fsmGuiScanDataBase(fsmExpParamPath);
    
    % Fill the scroll-down menu in the user interface
    labels=cell(length(fsmExpParam)+1,1);
    labels(1,1)={'Select experiment'};
    for i=1:length(fsmExpParam)
        labels(i+1,1)={fsmExpParam(i).label};
    end
    set(handles.expPopup,'String',labels);
    
    % Attach fsmExpParam to the userData field of the scroll-down menu
    set(handles.expPopup,'UserData',fsmExpParam);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% START CALCULATION CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = start_Callback(h, eventdata, handles, varargin)

% Read fsmParam from start button
fsmParam=get(handles.start,'UserData');

% Complement fsmParam with values from GUI
[fsmParam,status]=fsmGuiReadParameters(handles,fsmParam);
if status==0
    return
end

% Keep children status
allChildren=allchild(handles.fsmGuiMain);
childrenStatus=get(allChildren,'Enable');

% Disable settings window
set(allChildren,'Enable','Off');

% Start main program
fsmParam=fsmMain(fsmParam);

% Update GUI
fsmGuiWriteParameters(fsmParam,handles);

% Read parameter structure from handles.expPopup
fsmExpParam=get(handles.expPopup,'UserData');

% Save selected settings to disk (file: settings.txt)
fsmWriteParamsToTextFile(fsmParam,fsmExpParam,fsmParam.main.path)

% Link modified fsmParam to start button
set(handles.start,'UserData',fsmParam);

% Re-enable settings window
for i=1:length(allChildren)
	set(allChildren(i),'Enable',char(childrenStatus(i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFAULT PARAMETERS CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = defaultButton_Callback(h, eventdata, handles, varargin)

% Read defaultFsmParam from default button
defaultFsmParam=get(handles.defaultButton,'UserData');

% Set default values
fsmGuiWriteParameters(defaultFsmParam,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONFIDENCE INTERVAL CALLBACKS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = confOne_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(1)));
set(handles.confOne,'Value',1);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = confTwo_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(2)));
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',1);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = confThree_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(3)));
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',1);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = confFour_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(4)));
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',1);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = confFive_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(5)));
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',1);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = confSix_Callback(h, eventdata, handles, varargin)

% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
set(handles.editZValue,'String',num2str(zValue(6)));
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BLEACHING CALLBACKS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = bleachRadioOff_Callback(h, eventdata, handles, varargin)

set(handles.bleachRadioOff,'Value',1);
set(handles.bleachRadio1x,'Value',0);
set(handles.bleachRadio2x,'Value',0);
set(handles.bleachRadio3x,'Value',0);

% --------------------------------------------------------------------
function varargout = bleachRadio1x_Callback(h, eventdata, handles, varargin)

set(handles.bleachRadioOff,'Value',0);
set(handles.bleachRadio1x,'Value',1);
set(handles.bleachRadio2x,'Value',0);
set(handles.bleachRadio3x,'Value',0);

% --------------------------------------------------------------------
function varargout = bleachRadio2x_Callback(h, eventdata, handles, varargin)

set(handles.bleachRadioOff,'Value',0);
set(handles.bleachRadio1x,'Value',0);
set(handles.bleachRadio2x,'Value',1);
set(handles.bleachRadio3x,'Value',0);

% --------------------------------------------------------------------
function varargout = bleachRadio3x_Callback(h, eventdata, handles, varargin)

set(handles.bleachRadioOff,'Value',0);
set(handles.bleachRadio1x,'Value',0);
set(handles.bleachRadio2x,'Value',0);
set(handles.bleachRadio3x,'Value',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRIANGULATION CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = TriangCheck_Callback(h, eventdata, handles, varargin)
% Nothing must be done here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTO POLYGON CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = autoPolCheck_Callback(h, eventdata, handles, varargin)
% Nothing has to be done

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PATH SELECTION CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = ButtonBrowse_Callback(h, eventdata, handles, varargin)

% Select path
matlabVersion=ver('MATLAB');
if str2num(matlabVersion.Version)<6.5
    [newfile,newpath]=uiputfile('fsm.ini','Select work directory');
else
    newpath=uigetdir('','Select work directory');
end

if newpath~=0
    set(handles.pathEdit,'String',newpath);

    % If an fsmParam.mat file is found in this path, use this as fsmParam
    if exist([newpath,filesep,'fsmParam.mat'])==2
        
        % Inform user
        uiwait(msgbox('A file ''fsmParam.mat'' has been found in the selected path. Loading and setting stored values.','Info','help','modal'));
        
        % Load the file found
        load([newpath,filesep,'fsmParam.mat']);  
            
        % Fill all fields with the values from fsmParam
        fsmGuiWriteParameters(fsmParam,handles);
        
        % Link fsmParam to the start button (the default value are not modified)
        set(handles.start,'UserData',fsmParam);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREPROCESSING MODULE ON/OFF CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = checkPrepModule_Callback(h, eventdata, handles, varargin)

if get(handles.checkPrepModule,'Value')==0
    set(handles.TriangCheck,'Enable','off');
    set(handles.autoPolCheck,'Enable','off');
    set(handles.textDel,'Enable','off');
    set(handles.textCameraCalPar,'Enable','off');
    set(handles.expPopup,'Enable','off');
    set(handles.textAdvancedPrep,'Enable','off');
    set(handles.textGauss,'Enable','off');
    set(handles.editGauss,'Enable','off');
    set(handles.textExplanation,'Enable','off');
    set(handles.textExplanation2,'Enable','off');   
    set(handles.primaryRadio,'Enable','off');
    set(handles.tertiaryRadio,'Enable','off');
else
    set(handles.TriangCheck,'Enable','on');
    set(handles.autoPolCheck,'Enable','on');
    set(handles.textDel,'Enable','on');
    set(handles.bleachRadioOff,'Enable','on');
    set(handles.textCameraCalPar,'Enable','on');
    set(handles.expPopup,'Enable','on');
    set(handles.expPopup,'Enable','on');
    set(handles.textAdvancedPrep,'Enable','on');
    set(handles.textGauss,'Enable','on');
    set(handles.editGauss,'Enable','on');
    set(handles.textExplanation,'Enable','on');
    set(handles.textExplanation2,'Enable','on');   
    set(handles.primaryRadio,'Enable','on');
    set(handles.tertiaryRadio,'Enable','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% KINETIC ANALYSIS MODULE ON/OFF CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = checkKinModule_Callback(h, eventdata, handles, varargin)

if get(handles.checkKinModule,'Value')==0
    set(handles.textBleach,'Enable','off');
    set(handles.bleachRadioOff,'Enable','off');
    set(handles.bleachRadio1x,'Enable','off');
    set(handles.bleachRadio2x,'Enable','off');
    set(handles.bleachRadio3x,'Enable','off');
else
    set(handles.textBleach,'Enable','on');
    set(handles.bleachRadioOff,'Enable','on');
    set(handles.bleachRadio1x,'Enable','on');
    set(handles.bleachRadio2x,'Enable','on');
    set(handles.bleachRadio3x,'Enable','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRACKER THRESHOLD CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = editThreshold_Callback(h, eventdata, handles, varargin)
% Nothing to do here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRACKER MODULE ON/OFF CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkTrackModule_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkTrackModule
if get(hObject,'Value')==0
    set(handles.radioTrackBrownian,'Enable','off');
    set(handles.radioEnhTrackBrownian,'Enable','off');
    set(handles.radioTrackFlow,'Enable','off');
    set(handles.textThreshold,'Enable','off');
    set(handles.editThreshold,'Enable','off');
    set(handles.checkEnhTrack,'Enable','off');
    set(handles.checkGrid,'Enable','off');
else
    set(handles.radioTrackBrownian,'Enable','on');
    set(handles.radioEnhTrackBrownian,'Enable','on');
    set(handles.radioTrackFlow,'Enable','on');
    set(handles.checkEnhTrack,'Enable','on');
    set(handles.textThreshold,'Enable','on');
    set(handles.editThreshold,'Enable','on');
    if get(handles.checkEnhTrack,'Value')==1
        set(handles.checkGrid,'Enable','on');
    else
        set(handles.checkGrid,'Enable','off');
    end               
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRACKER SELECTION CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radioTrackBrownian_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of radioTrackBrownian
set(handles.radioTrackBrownian,'Value',1);
set(handles.radioEnhTrackBrownian,'Value',0);
set(handles.radioTrackFlow,'Value',0);
% set(handles.textThreshold,'Enable','off');
% set(handles.editThreshold,'Enable','off');

function radioEnhTrackBrownian_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radioEnhTrackBrownian
set(handles.radioTrackBrownian,'Value',0);
set(handles.radioEnhTrackBrownian,'Value',1);
set(handles.radioTrackFlow,'Value',0);
% set(handles.textThreshold,'Enable','on');
% set(handles.editThreshold,'Enable','on');

% -------------------------------------------------------------------------

function radioTrackFlow_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radioTrackFlow
set(handles.radioTrackBrownian,'Value',0);
set(handles.radioEnhTrackBrownian,'Value',0);
set(handles.radioTrackFlow,'Value',1);
set(handles.textThreshold,'Enable','on');
set(handles.editThreshold,'Enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SAVE PARAMETERS (FOR BATCH JOBS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveParamsButton_Callback(hObject, eventdata, handles)
userPath=get(handles.pathEdit,'String');
if isempty(userPath)
    uiwait(msgbox('Please specify a work path.','Warning','warn'));
    return;
else
    % Create ROOT directory
    if exist(userPath)~=7
        % Directory does not exist - create it
        if userPath(2)==':';   % Windows
            % Drive letter specified
            st=mkdir(userPath(1:3),userPath(4:end));
        else
            st=mkdir(userPath);
        end
    end

    % Read fsmParam from start button
    fsmParam=get(handles.start,'UserData');
    % Read parameters from GUI
    [fsmParam,status]=fsmGuiReadParameters(handles,fsmParam);
    if status==0
        return
    end

    % Save parameters
    olddir=cd;
    cd(userPath);
    save fsmParam.mat fsmParam;
    cd(olddir);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RESULT DISPLAY CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in checkDispModule.
function checkDispModule_Callback(hObject, eventdata, handles)
% Nothing to do

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BUILDER CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in checkBuildModule.
function checkBuildModule_Callback(hObject, eventdata, handles)
% Nothing to do


% --- Executes during object creation, after setting all properties.
function editGass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editGass_Callback(hObject, eventdata, handles)
% hObject    handle to editGass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGass as text
%        str2double(get(hObject,'String')) returns contents of editGass as a double


% --- Executes during object creation, after setting all properties.
function editZValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editZValue_Callback(hObject, eventdata, handles)

% If  a z-value is entered by the user, turn all radio button to 0
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% --------------------------------------------------------------------
function varargout = primaryRadio_Callback(h, eventdata, handles, varargin)

set(handles.primaryRadio,'Value',1);
set(handles.tertiaryRadio,'Value',0);
set(handles.scaleRadio,'Value',0);
set(handles.percText,'Enable','off');
set(handles.percEdit,'Enable','off');
set(handles.orderText,'Enable','off');
set(handles.orderEdit,'Enable','off');
set(handles.sigText,'Enable','off');
set(handles.sigEdit,'Enable','off');
% Turn on kinetic analysis
kineticAnalysis(handles,'on');

% --------------------------------------------------------------------
function varargout = tertiaryRadio_Callback(h, eventdata, handles, varargin)

set(handles.primaryRadio,'Value',0);
set(handles.tertiaryRadio,'Value',1);
set(handles.scaleRadio,'Value',0);
set(handles.percText,'Enable','on');
set(handles.percEdit,'Enable','on');
set(handles.orderText,'Enable','on');
set(handles.orderEdit,'Enable','on');
set(handles.sigText,'Enable','off');
set(handles.sigEdit,'Enable','off');
% Turn on kinetic analysis
kineticAnalysis(handles,'on');

% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on button press in checkEnhTrack.
function checkEnhTrack_Callback(hObject, eventdata, handles)
% hObject    handle to checkEnhTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkEnhTrack
if get(handles.checkEnhTrack,'Value')==1
    set(handles.checkGrid,'Enable','on');
else
    set(handles.checkGrid,'Enable','off');
end    

% --- Executes on selection change in expPopup.
function expPopup_Callback(hObject, eventdata, handles)
% hObject    handle to expPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns expPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from expPopup

% Read experiment value
exp=get(handles.expPopup,'Value');

% Read parameter structure from handles.expPopup
fsmExpParam=get(handles.expPopup,'UserData');

% Read defaultFsmParam from 'UserData' field of default button
defaultFsmParam=get(handles.defaultButton,'UserData');

% Write bit depth
if exp==1
    set(handles.bitDepthEdit,'String',num2str(log2(defaultFsmParam.main.normMax+1)));    
else
    set(handles.bitDepthEdit,'String',num2str(fsmExpParam(exp-1).bitDepth));
end

% Write gauss ratio
if exp==1
    set(handles.editGauss,'String',num2str(defaultFsmParam.prep.gaussRatio));
else
    set(handles.editGauss,'String',num2str(fsmExpParam(exp-1).gaussRatio));
end

% Write description
if exp==1
    set(handles.textDescr,'String','Experiment description');
else
    set(handles.textDescr,'String',fsmExpParam(exp-1).description);
end


function orderEdit_Callback(hObject, eventdata, handles)
% hObject    handle to orderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of orderEdit as text
%        str2double(get(hObject,'String')) returns contents of orderEdit as a double

function percEdit_Callback(hObject, eventdata, handles)
% hObject    handle to percEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of percEdit as text
%        str2double(get(hObject,'String')) returns contents of percEdit as a double


% --- Executes during object creation, after setting all properties.
function percEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function orderEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --------------------------------------------------------------------
function varargout = sigText_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = sigEdit_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = scaleRadio_Callback(h, eventdata, handles, varargin)

set(handles.primaryRadio,'Value',0);
set(handles.tertiaryRadio,'Value',0);
set(handles.scaleRadio,'Value',1);
set(handles.percText,'Enable','off');
set(handles.percEdit,'Enable','off');
set(handles.orderText,'Enable','off');
set(handles.orderEdit,'Enable','off');
set(handles.sigText,'Enable','on');
set(handles.sigEdit,'Enable','on');
% Turn off kinetic analysis
kineticAnalysis(handles,'off');

% --------------------------------------------------------------------
function kineticAnalysis(handles,flag)

if ~strcmp(flag,'off') & ~strcmp(flag,'on')
    error('Wrong value for parameter flag');
end

% Enable/Disable
set(handles.checkBuildModule,'Enable',flag); % Builder module
set(handles.checkKinModule,'Enable',flag);   % Kinetic analysis module
set(handles.textBleach,'Enable',flag);
set(handles.bleachRadioOff,'Enable',flag);
set(handles.bleachRadio1x,'Enable',flag);
set(handles.bleachRadio2x,'Enable',flag);
set(handles.bleachRadio3x,'Enable',flag);
set(handles.checkDispModule,'Enable',flag);  % Result display module

% If flag is 'off', it means that the scale space approach has been selected.
% If this is the case, not only must the build, kin, and disp module be disabled, such
% that they cannot be selected by the user, but they also have to be turned off such 
% that the software is not going to run them. The opposite is not necessary

if strcmp(flag,'off')
    
    set(handles.checkBuildModule,'Value',0); % Builder module
    set(handles.checkKinModule,'Value',0);   % Kinetic analysis module
    set(handles.checkDispModule,'Value',0);  % Result display module
    
end


% --- Executes on button press in checkGrid.
function checkGrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkGrid
if get(handles.checkGrid,'Value')==1
    uiwait(msgbox('Use this field only in case of OUT OF MEMORY problems during tracking.','Not recommended','help','modal'));
end

% --- Executes on button press in radioDispModeOrig.
function radioDispModeOrig_Callback(hObject, eventdata, handles)
% hObject    handle to radioDispModeOrig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioDispModeOrig
set(handles.radioDispModeOrig,'Value',1); % Keep it on for the moment


% --- Executes on button press in radioDispModeTCO.
function radioDispModeTCO_Callback(hObject, eventdata, handles)
% hObject    handle to radioDispModeTCO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioDispModeTCO


% --- Executes on button press in drawROICheck.
function drawROICheck_Callback(hObject, eventdata, handles)
% hObject    handle to drawROICheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawROICheck


% --------------------------------------------------------------------
function fsmGuiMain_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fsmGMH=findall(0,'Tag','fsmGuiMain'); % Get the handle of fsmGuiMain
choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
switch choice,
    case 'Yes', delete(fsmGMH);
    case 'No', return;
end % switch


% --- Executes during object deletion, before destroying properties.
function fsmGuiMain_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fsmGuiMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


