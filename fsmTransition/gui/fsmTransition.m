function varargout = fsmTransition(varargin)
% FSMTRANSITION M-file for fsmTransition.fig
%      FSMTRANSITION, by itself, creates a new FSMTRANSITION or raises the existing
%      singleton*.
%
%      H = FSMTRANSITION returns the handle to a new FSMTRANSITION or the handle to
%      the existing singleton*.
%
%      FSMTRANSITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMTRANSITION.M with the given input arguments.
%
%      FSMTRANSITION('Property','Value',...) creates a new FSMTRANSITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmTransition_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmTransition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmTransition

% Last Modified by GUIDE v2.5 03-Sep-2004 14:57:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmTransition_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmTransition_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fsmTransition is made visible.
function fsmTransition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmTransition (see VARARGIN)

% Choose default command line output for fsmTransition
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Update project info from fsmCenter
updateProjectInfo(handles);


% UIWAIT makes fsmTransition wait for user response (see UIRESUME)
% uiwait(handles.fsmTransition);


% --- Outputs from this function are returned to the command line.
function varargout = fsmTransition_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CREATEFCN
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function popupTack_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editPanelProfiles_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  READ PRPANEL FILES
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkPanelProfiles_Callback(hObject, eventdata, handles)

switch get(handles.checkPanelProfiles,'Value')
    case 0, set(handles.editPanelProfiles,'Enable','off'); set(handles.textPanel1,'Enable','off');
    case 1, set(handles.editPanelProfiles,'Enable','on'); set(handles.textPanel1,'Enable','on');
    otherwise, error('Field editPanelProfiles has value different from 0 | 1');
end

function checkPanelNormals_Callback(hObject, eventdata, handles)
switch get(handles.checkPanelNormals,'Value')
    case 0, set(handles.checkPanelProfiles,'Enable','off');set(handles.editPanelProfiles,'Enable','off'); set(handles.textPanel1,'Enable','off');
    case 1, set(handles.checkPanelProfiles,'Enable','on');set(handles.editPanelProfiles,'Enable','on'); set(handles.textPanel1,'Enable','on');
    otherwise, error('Field editPanelProfiles has value different from 0 | 1');
end

function checkPanelProtrusions_Callback(hObject, eventdata, handles)

function checkPanelEdgePixels_Callback(hObject, eventdata, handles)

function editPanelProfiles_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMPORT PRPANEL DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushPanel_Callback(hObject, eventdata, handles)
% Get project dir
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal');
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    errordlg('Could not find all (sub)projects. Quitting.','Error','modal');
    return
end

% Form path for the protrusion subdirectory
protrusionPath=[edgeProjDir,filesep,'protrusion'];

% Check which of the checkboxes are on
switch get(handles.checkPanelEdgePixels,'Value') % Edge pixels
    case 0, edgePixelsFileName=-1;
    case 1, edgePixelsFileName=[protrusionPath,filesep,'pixel_edge.dat'];
    otherwise, error('Field checkPanelEdgePixels has a wrong value');
end
switch get(handles.checkPanelNormals,'Value') % Normals
    case 0, normalsFileName=-1;
    case 1, normalsFileName=[protrusionPath,filesep,'av_normal.dat'];
    otherwise, error('Field checkPanelNormals has a wrong value');
end
switch get(handles.checkPanelProfiles,'Value') % Profiles
    case 0, toggleProfiles=0;
    case 1, toggleProfiles=1;
    otherwise, error('Field checkPanelProfiles has a wrong value');
end
switch get(handles.checkPanelProtrusions,'Value') % Protrusions
    case 0, protFileName=-1;
    case 1, protFileName=[protrusionPath,filesep,'av_prot.dat'];
    otherwise, error('Field checkPanelProtrusions has a wrong value');
end
vLength=str2num(get(handles.editPanelProfiles,'String')); % Profile lengths
if normalsFileName==-1
    toggleProfiles=0; % If normals is turned off, no profiles can be calculated
end
switch get(handles.checkPanelSegments,'Value') % Segments
    case 0, protFileName=-1;
    case 1, 
        d=dir(protrusionPath);
        segmFileName=[];
        for i=1:length(d)
            if strncmp(d(i).name,'s_mask_',7)==1
                segmFileName=[protrusionPath,filesep,d(i).name];
                break;
            end
        end
    otherwise, error('Field checkPanelSegments has a wrong value');
end
if isempty(segmFileName)
    disp('Files for segments could not be found.');
end
% Run function
[normals,profiles,protrusions,edgePixels,segments]=fsmTransReadPrPanelDatFiles(edgePixelsFileName,normalsFileName,toggleProfiles,protFileName,segmFileName,vLength);
% Save outputs to disk
save([lplaProjDir,filesep,'normals.mat'],'normals');
save([lplaProjDir,filesep,'profiles.mat'],'profiles');
save([lplaProjDir,filesep,'protrusions.mat'],'protrusions');
save([lplaProjDir,filesep,'edgePixels.mat'],'edgePixels');
save([lplaProjDir,filesep,'segments.mat'],'segments');
msg=['Files saved in [',lplaProjDir,']'];
uiwait(msgbox(msg,'Info','modal'));
% Return outputs to Matlab workspace
assignin('base','normals',normals);
assignin('base','profiles',profiles);
assignin('base','protrusions',protrusions);
assignin('base','edgePixels',edgePixels);
assignin('base','segments',segments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END IMPORT PRPANEL DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  EXTRACT TRANSITION
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushExtract_Callback(hObject, eventdata, handles)
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal');
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    errordlg('Could not find all (sub)projects. Quitting.','Error','modal');
    return
end

toggleKinetic=get(handles.checkExtractKin,'Value');
toggleSpeed=get(handles.checkExtractSpeed,'Value');

% Load profiles.mat
profilesFileName=[lplaProjDir,filesep,'profiles.mat'];
if exist(profilesFileName,'file')==0
    errordlg('Could not find file ''profiles.mat'' in your lpla subdirectory. Make sure to ''Import prPanel data''.','Error','modal');
    return
end
try
    s=load(profilesFileName);
    profiles=s.profiles;
catch
    errordlg('Corrupted ''profiles.mat'' file. Try to re-import it (''Import prPanel data'') or re-run prPanel .','Error','modal');
    return
end

n=str2num(get(handles.editExtractFrames,'String'));
medianFiltFrames=str2num(get(handles.editExtractMedian,'String'));
dist=str2num(get(handles.editExtractSupport,'String'));
toggleDpRatio=get(handles.checkExtractDpRatio,'Value');
[positionsScore,positionsVel,positionsTrans,fPositionsTrans,dpRatio,scoreProfiles,speedProfiles]=fsmTransExtractLpToLaTransition(tackProjDir,toggleKinetic,toggleSpeed,toggleDpRatio,profiles,dist,n,medianFiltFrames);

eval(['save ',lplaProjDir,filesep,'positionsScore.mat positionsScore']);
eval(['save ',lplaProjDir,filesep,'positionsVel.mat positionsVel']);
eval(['save ',lplaProjDir,filesep,'positionsTrans.mat positionsTrans']);
eval(['save ',lplaProjDir,filesep,'fPositionsTrans.mat fPositionsTrans']);
eval(['save ',lplaProjDir,filesep,'dpRatio.mat dpRatio']);
eval(['save ',lplaProjDir,filesep,'scoreProfiles.mat scoreProfiles']);
eval(['save ',lplaProjDir,filesep,'speedProfiles.mat speedProfiles']);
% Save settings
fid=fopen([lplaProjDir,filesep,'extractTransitionSettings.txt'],'a+');
if fid~=-1
    fprintf(fid,'>%s\n\n',datestr(now,0));
    fprintf(fid,'Use kinetic criterion      : %d\n',toggleKinetic);
    fprintf(fid,'Use speed criterion        : %d\n',toggleSpeed);
    fprintf(fid,'Calculate depoly/poly ratio: %d\n',toggleDpRatio);
    fprintf(fid,'Median filter order        : %d\n',medianFiltFrames);
    fprintf(fid,'Support (pixels)           : %d\n',dist);
    fprintf(fid,'Time averaging (frames)    : %d\n\n',n);
    fclose(fid);
end
msg=['Files saved in [',lplaProjDir,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','positionsScore',positionsScore);
assignin('base','positionsVel',positionsVel);
assignin('base','positionsTrans',positionsTrans);
assignin('base','fPositionsTrans',fPositionsTrans);
assignin('base','dpRatio',dpRatio);
assignin('base','scoreProfiles',scoreProfiles);
assignin('base','speedProfiles',speedProfiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END EXTRACT TRANSITION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in checkExtractKin.
function checkExtractKin_Callback(hObject, eventdata, handles)
% hObject    handle to checkExtractKin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkExtractKin
if get(handles.checkExtractKin,'Value')==1
    set(handles.checkExtractDpRatio,'Enable','on');
else
    set(handles.checkExtractDpRatio,'Enable','off');
end    

% --- Executes on button press in checkExtractSpeed.
function checkExtractSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to checkExtractSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkExtractSpeed


% --- Executes during object creation, after setting all properties.
function editExtractMedian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editExtractMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editExtractMedian_Callback(hObject, eventdata, handles)
% hObject    handle to editExtractMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editExtractMedian as text
%        str2double(get(hObject,'String')) returns contents of editExtractMedian as a double


% --- Executes during object creation, after setting all properties.
function editExtractFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editExtractFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editExtractFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editExtractFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editExtractFrames as text
%        str2double(get(hObject,'String')) returns contents of editExtractFrames as a double


% --- Executes on button press in pushLpPolygons.
function pushLpPolygons_Callback(hObject, eventdata, handles)
% hObject    handle to pushLpPolygons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
method=find([get(handles.radioLpMethod1,'Value') get(handles.radioLpMethod2,'Value')]);
if method==1
    % Select splineCoordsTrans.mat
    [fName,dirName] = uigetfile(...
        {'splineCoordsTrans.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select splineCoordsTrans.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    coordsTrans=splineCoordsTrans;
    % Select splineCoordsEdge.mat
    [fName,dirName] = uigetfile(...
        {'splineCoordsEdge.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select splineCoordsEdge.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    coordsEdge=splineCoordsEdge;
else
    % Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat
    [fName,dirName] = uigetfile(...
        {'*.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    trans=load([dirName,fName]);
    field=fieldnames(trans);
    if strcmp(field,'fPositionsTrans')==0 & strcmp(field,'positionsTrans')==0 & strcmp(field,'positionsScore')==0 & strcmp(field,'positionsVel')==0 
        error('The loaded file must be one of fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
    end
    eval(['coordsTrans=trans.',char(field),';']);
    % Select edgePixels.mat
    [fName,dirName] = uigetfile(...
        {'edgePixels.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select edgePixels.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    coordsEdge=edgePixels;
end
% Frames used to extract transition
frames=str2num(get(handles.editLpFrames,'String'));
lpPolygons=createLpPolygons(method,coordsEdge,coordsTrans,frames,get(handles.checkLpShowPolygons,'Value'));
% Save output
path=uigetdir('','Select directory where to save the computed polygons');
if path~=0
    eval(['save ',path,filesep,'lpPolygons.mat lpPolygons']);
    msg=['Files saved in [',path,']'];
    uiwait(msgbox(msg,'Info','modal'));
end
% Return to Matlab workspace
assignin('base','lpPolygons',lpPolygons);


% --- Executes on button press in radioLpMethod1.
function radioLpMethod1_Callback(hObject, eventdata, handles)
% hObject    handle to radioLpMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioLpMethod1
set(handles.radioLpMethod1,'Value',1);
set(handles.radioLpMethod2,'Value',0);

% --- Executes on button press in radioLpMethod2.
function radioLpMethod2_Callback(hObject, eventdata, handles)
% hObject    handle to radioLpMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioLpMethod2
set(handles.radioLpMethod2,'Value',1);
set(handles.radioLpMethod1,'Value',0);


% --- Executes on button press in checkLpShowPolygons.
function checkLpShowPolygons_Callback(hObject, eventdata, handles)
% hObject    handle to checkLpShowPolygons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkLpShowPolygons


% --- Executes on button press in pushTransToEdgeDist.
function pushTransToEdgeDist_Callback(hObject, eventdata, handles)
% hObject    handle to pushTransToEdgeDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat
[fName,dirName] = uigetfile(...
    {'*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
trans=load([dirName,fName]);
field=fieldnames(trans);
if strcmp(field,'fPositionsTrans')==0 & strcmp(field,'positionsTrans')==0 & strcmp(field,'positionsScore')==0 & strcmp(field,'positionsVel')==0 
    error('The loaded file must be one of fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
end
eval(['positionsTrans=trans.',char(field),';']);
% Select edgePixels.mat
[fName,dirName] = uigetfile(...
    {'edgePixels.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select edgePixels.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
load([dirName,fName]);
splineTolEdge=str2num(get(handles.editTEDTolEdge,'String'));
splineTolTrans=str2num(get(handles.editTEDTolTrans,'String'));
medianFilterFrames=str2num(get(handles.editTEDTolMedian,'String'));
[transToEdgeDistNearest,transToEdgeDistMech,spTrans_y,spTrans_x,spEdge_y,spEdge_x,splineCoordsTrans,splineCoordsEdge]=fsmTransTransitionToEdgeDistances(positionsTrans,edgePixels,splineTolEdge,splineTolTrans,medianFilterFrames);
% Save output
path=uigetdir('','Select directory where to save the computed distances');
if path~=0
    eval(['save ',path,filesep,'transToEdgeDistNearest.mat transToEdgeDistNearest']);
    eval(['save ',path,filesep,'transToEdgeDistMech.mat transToEdgeDistMech']);
    eval(['save ',path,filesep,'spEdge.mat spEdge_y spEdge_x']);
    eval(['save ',path,filesep,'spTrans.mat spTrans_y spTrans_x']);
    eval(['save ',path,filesep,'splineCoordsTrans.mat splineCoordsTrans']);
    eval(['save ',path,filesep,'splineCoordsEdge.mat splineCoordsEdge']);
    eval(['save ',path,filesep,'settings.mat splineTolEdge splineTolTrans medianFilterFrames']);
    msg=['Files saved in [',path,']'];
    uiwait(msgbox(msg,'Info','modal'));
end


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function editExtractSupport_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editExtractSupport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editExtractSupport_Callback(hObject, eventdata, handles)
% hObject    handle to editExtractSupport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editExtractSupport as text
%        str2double(get(hObject,'String')) returns contents of editExtractSupport as a double


% --- Executes on button press in checkExtractdpRatio.
function checkExtractdpRatio_Callback(hObject, eventdata, handles)
% hObject    handle to checkExtractdpRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkExtractdpRatio


% --- Executes on button press in checkExtractDpRatio.
function checkExtractDpRatio_Callback(hObject, eventdata, handles)
% hObject    handle to checkExtractDpRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkExtractDpRatio


% --- Executes during object creation, after setting all properties.
function editLpFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLpFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editLpFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editLpFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLpFrames as text
%        str2double(get(hObject,'String')) returns contents of editLpFrames as a double


% --- Executes on button press in pushPlotLpPolygons.
function pushPlotLpPolygons_Callback(hObject, eventdata, handles)
% hObject    handle to pushPlotLpPolygons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load lpPolygons
[fName,dirName] = uigetfile(...
    {'lpPolygons.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select lpPolygons.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
lp=load([dirName,fName]);
field=fieldnames(lp);
eval(['lpPolygons=lp.',char(field),';']);
plotLpPolygons(lpPolygons);


% --- Executes on button press in pushSpline.
function pushSpline_Callback(hObject, eventdata, handles)
% hObject    handle to pushSpline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat
[fName,dirName] = uigetfile(...
    {'*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
trans=load([dirName,fName]);
field=fieldnames(trans);
if strcmp(field,'fPositionsTrans')==0 & strcmp(field,'positionsTrans')==0 & strcmp(field,'positionsScore')==0 & strcmp(field,'positionsVel')==0 
    error('The loaded file must be one of fPositionsScore.mat or positionsScore.mat or positionsScore.mat or positionsVel.mat');
end
eval(['positionsTrans=trans.',char(field),';']);
% Select edgePixels.mat
[fName,dirName] = uigetfile(...
    {'edgePixels.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select edgePixels.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
load([dirName,fName]);
toleranceEdge=str2num(get(handles.editSplineTolEdge,'String'));
toleranceTrans=str2num(get(handles.editSplineTolTrans,'String'));
medianFilterFrames=str2num(get(handles.editSplineMedian,'String'));
points=str2num(get(handles.editSplinePoints,'String'));
[spEdge_y,spEdge_x,spTrans_y,spTrans_x,splineCoordsEdge,splineCoordsTrans]=fsmTransFitSplines(positionsTrans,edgePixels,medianFilterFrames,toleranceTrans,toleranceEdge,points)
% Save output
path=uigetdir('','Select directory where to save the computed splines');
if path~=0
    eval(['save ',path,filesep,'spEdge.mat spEdge_y spEdge_x']);
    eval(['save ',path,filesep,'spTrans.mat spTrans_y spTrans_x']);
    eval(['save ',path,filesep,'splineCoordsEdge.mat splineCoordsEdge']);
    eval(['save ',path,filesep,'splineCoordsTrans.mat splineCoordsTrans']);
    eval(['save ',path,filesep,'settings.mat toleranceEdge toleranceTrans medianFilterFrames']);
    msg=['Files saved in [',path,']'];
    uiwait(msgbox(msg,'Info','modal'));
end


% --- Executes during object creation, after setting all properties.
function editSplineTolEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSplineTolEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSplineTolEdge_Callback(hObject, eventdata, handles)
% hObject    handle to editSplineTolEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSplineTolEdge as text
%        str2double(get(hObject,'String')) returns contents of editSplineTolEdge as a double


% --- Executes during object creation, after setting all properties.
function editSplineTolTrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSplineTolTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSplineTolTrans_Callback(hObject, eventdata, handles)
% hObject    handle to editSplineTolTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSplineTolTrans as text
%        str2double(get(hObject,'String')) returns contents of editSplineTolTrans as a double



% --- Executes during object creation, after setting all properties.
function editSplineMedian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSplineMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSplineMedian_Callback(hObject, eventdata, handles)
% hObject    handle to editSplineMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSplineMedian as text
%        str2double(get(hObject,'String')) returns contents of editSplineMedian as a double


% --- Executes on button press in pushDt.
function pushDt_Callback(hObject, eventdata, handles)
% hObject    handle to pushDt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fName,dirName] = uigetfile(...
    {'spEdge.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select spEdge.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
load([dirName,fName]);

[fName,dirName] = uigetfile(...
    {'spTrans.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select spTrans.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
load([dirName,fName]);

points=1:str2num(get(handles.editDtPoints,'String'));

[distancesNearestEdge,distancesMechEdge]=fsmTransDistancesDt(spEdge_y,spEdge_x,points);
[distancesNearestTrans,distancesMechTrans]=fsmTransDistancesDt(spTrans_y,spTrans_x,points);

% Save output
path=uigetdir('','Select directory where to save the computed protrusions');
if path~=0
    eval(['save ',path,filesep,'distancesNearestEdge.mat distancesNearestEdge']);
    eval(['save ',path,filesep,'distancesMechEdge.mat distancesMechEdge']);
    eval(['save ',path,filesep,'distancesNearestTrans.mat distancesNearestTrans']);
    eval(['save ',path,filesep,'distancesMechTrans.mat distancesMechTrans']);
    eval(['save ',path,filesep,'points.mat points']);
    msg=['Files saved in [',path,']'];
    uiwait(msgbox(msg,'Info','modal'));
end
% Return to Matlab workspace
assignin('base','distancesNearestEdge',distancesNearestEdge);
assignin('base','distancesMechEdge',distancesMechEdge);
assignin('base','distancesNearestTrans',distancesNearestTrans);
assignin('base','distancesMechTrans',distancesMechTrans);
assignin('base','points',points);

% --- Executes during object creation, after setting all properties.
function editDtPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDtPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editDtPoints_Callback(hObject, eventdata, handles)
% hObject    handle to editDtPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDtPoints as text
%        str2double(get(hObject,'String')) returns contents of editDtPoints as a double


% --- Executes on button press in pushProtrusionMaps.
function pushProtrusionMaps_Callback(hObject, eventdata, handles)
% hObject    handle to pushProtrusionMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fName,dirName] = uigetfile(...
    {'transToEdgeDist*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select transToEdgeDistNearest.mat or transToEdgeDistMech.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
trans=load([dirName,fName]);
field=fieldnames(trans);
if strcmp(field,'transToEdgeDistNearest')==0 & strcmp(field,'transToEdgeDistMech')==0
    error('The loaded file must be one of transToEdgeDistNearest.mat or transToEdgeDistMech.mat');
end
eval(['transToEdgeDist=trans.',char(field),';']);
% Select distancesNearestEdge.mat or distancesMechEdge.mat
[fName,dirName] = uigetfile(...
    {'*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select distancesNearestEdge.mat or distancesMechEdge.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
trans=load([dirName,fName]);
field=fieldnames(trans);
if strcmp(field,'distancesNearestEdge')==0 & strcmp(field,'distancesMechEdge')==0
    error('The loaded file must be one of distancesNearestEdge.mat or distancesMechEdge.mat');
end
eval(['distancesEdge=trans.',char(field),';']);
% Select distancesNearestTrans.mat or distancesMechTrans.mat
[fName,dirName] = uigetfile(...
    {'*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select distancesNearestTrans.mat or distancesMechTrans.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
trans=load([dirName,fName]);
field=fieldnames(trans);
if strcmp(field,'distancesNearestTrans')==0 & strcmp(field,'distancesMechTrans')==0
    error('The loaded file must be one of distancesNearestTrans.mat or distancesMechTrans.mat');
end
eval(['distancesTrans=trans.',char(field),';']);
time=str2num(get(handles.editProtrusionMapsFrames,'String'));
interpStepY=str2num(get(handles.editProtrusionMapsInterpY,'String'));
interpStepX=str2num(get(handles.editProtrusionMapsInterpX,'String'));
[transProtrusionMap2D,edgeProtrusionMap2D,stdTransProtrusionMap2D,stdEdgeProtrusionMap2D,indices]=fsmTransProtrusionMaps2d(transToEdgeDist,distancesTrans,distancesEdge,time,[interpStepY interpStepX]);
% Save output
path=uigetdir('','Select directory where to save the computed protrusion maps');
if path~=0
    eval(['save ',path,filesep,'transProtrusionMap2D.mat transProtrusionMap2D']);
    eval(['save ',path,filesep,'edgeProtrusionMap2D.mat edgeProtrusionMap2D']);
    eval(['save ',path,filesep,'stdTransProtrusionMap2D.mat stdTransProtrusionMap2D']);
    eval(['save ',path,filesep,'stdEdgeProtrusionMap2D.mat stdEdgeProtrusionMap2D']);
    eval(['save ',path,filesep,'indices.mat indices']);
    msg=['Files saved in [',path,']'];
    uiwait(msgbox(msg,'Info','modal'));
end

% Return to Matlab workspace
assignin('base','transProtrusionMap2D',transProtrusionMap2D);
assignin('base','edgeProtrusionMap2D',edgeProtrusionMap2D);
assignin('base','stdTransProtrusionMap2D',stdTransProtrusionMap2D);
assignin('base','stdEdgeProtrusionMap2D',stdEdgeProtrusionMap2D);
assignin('base','indices',indices);

% --- Executes during object creation, after setting all properties.
function editMapsFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMapsFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editMapsFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editMapsFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMapsFrames as text
%        str2double(get(hObject,'String')) returns contents of editMapsFrames as a double


% --- Executes during object creation, after setting all properties.
function editSplinePoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSplinePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSplinePoints_Callback(hObject, eventdata, handles)
% hObject    handle to editSplinePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSplinePoints as text
%        str2double(get(hObject,'String')) returns contents of editSplinePoints as a double


% --- Executes during object creation, after setting all properties.
function editTEDTolEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTEDTolEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editTEDTolEdge_Callback(hObject, eventdata, handles)
% hObject    handle to editTEDTolEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTEDTolEdge as text
%        str2double(get(hObject,'String')) returns contents of editTEDTolEdge as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function editTEDTolMedian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTEDTolMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editTEDTolMedian_Callback(hObject, eventdata, handles)
% hObject    handle to editTEDTolMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTEDTolMedian as text
%        str2double(get(hObject,'String')) returns contents of editTEDTolMedian as a double


% --- Executes during object creation, after setting all properties.
function editTEDTolTrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTEDTolTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editTEDTolTrans_Callback(hObject, eventdata, handles)
% hObject    handle to editTEDTolTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTEDTolTrans as text
%        str2double(get(hObject,'String')) returns contents of editTEDTolTrans as a double


% --- Executes during object creation, after setting all properties.
function editProtrusionMapsFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editProtrusionMapsFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editProtrusionMapsFrames as text
%        str2double(get(hObject,'String')) returns contents of editProtrusionMapsFrames as a double


% --- Executes during object creation, after setting all properties.
function editProtrusionMapsInterpX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsInterpX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editProtrusionMapsInterpX_Callback(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsInterpX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editProtrusionMapsInterpX as text
%        str2double(get(hObject,'String')) returns contents of editProtrusionMapsInterpX as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function editProtrusionMapsInterpY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsInterpY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editProtrusionMapsInterpY_Callback(hObject, eventdata, handles)
% hObject    handle to editProtrusionMapsInterpY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editProtrusionMapsInterpY as text
%        str2double(get(hObject,'String')) returns contents of editProtrusionMapsInterpY as a double


% --- Executes on button press in checkPanelSegments.
function checkPanelSegments_Callback(hObject, eventdata, handles)
% hObject    handle to checkPanelSegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkPanelSegments


% --- Executes during object creation, after setting all properties.
function popupEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupEdge.
function popupEdge_Callback(hObject, eventdata, handles)
% hObject    handle to popupEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupEdge contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupEdge




function popupTack_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ACCESSORY FUNCTIONS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles)
% Initialize output
projDir=[]; tackProjDir=[]; edgeProjDir=[]; lplaProjDir=[];
% Chech whether a project exists
projDir=get(handles.textCurrentProject,'String');
if isempty(projDir)
    return
else
    
    if ~isdir(projDir)
        warndlg('The directory specified does not exist. Please check your project.','Warning','modal');
        projDir=[]; % Return empty
    end

    % TACK
    
    % Check for multiple subprojects
    subProjs=get(handles.popupTack,'String');
    subProj=subProjs(get(handles.popupTack,'Value'));
    % Check that the project directory exists
    tackProjDir=[projDir,filesep,char(subProj),filesep];
    if ~isdir(tackProjDir)
        warndlg('The directory specified does not exist. Please check your project.','Warning','modal');
        tackProjDir=[]; % Return empty
    end
    
    % EDGE
    
    % Check for multiple subprojects
    subProjs=get(handles.popupEdge,'String');
    subProj=subProjs(get(handles.popupEdge,'Value'));
    % Check that the project directory exists
    edgeProjDir=[projDir,filesep,char(subProj),filesep];
    if ~isdir(edgeProjDir)
        warndlg('The directory specified does not exist. Please check your project.','Warning','modal');
        edgeProjDir=[]; % Return empty
    end

    % LPLA
    
    % Check for multiple subprojects
    subProjs=get(handles.popupLpla,'String');
    subProj=subProjs(get(handles.popupLpla,'Value'));
    % Check that the project directory exists
    lplaProjDir=[projDir,filesep,char(subProj),filesep];
    if ~isdir(lplaProjDir)
        warndlg('The directory specified does not exist. Please check your project.','Warning','modal');
        lplaProjDir=[]; % Return empty
    end

end

%%%

function updateProjectInfo(handles);

% Check whether fsmCenter is running - if not, open it
hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if isempty(hfsmC)
    hfsmC=fsmCenter;
end
% Get current project
handlesFsmCenter=guidata(hfsmC);
projDir=get(handlesFsmCenter.textCurrentProject,'String');
% Get settings from fsmCenter
settings=get(handlesFsmCenter.fsmCenter,'UserData');
% Check
if ~isempty(settings)
    if ~strcmp(projDir,settings.projDir)
        error('projDir stored in fsmCenter''s UserData and in fsmCenter''s .textCurrentProject field are different!');    
    end
end

set(handles.textCurrentProject,'String',projDir);

% If necessary, look for subproject
if isempty(projDir)
    set(handles.popupTack,'Enable','off');
    set(handles.popupEdge,'Enable','off');
else
    set(handles.textCurrentProject,'String',projDir);

    % Subproject tack
    subProjectsTack=findProjSubDir(projDir,'tack');
    set(handles.popupTack,'String',subProjectsTack);
    % Set as default the one chosen in the project
    [found,pos]=ismember(settings.subProjects(1),subProjectsTack);
    if found==0
        error('The directory stored in fsmCenter''s settings does not actually exist on disk!');
    else
        set(handles.popupTack,'Value',pos);
    end
    set(handles.popupTack,'Enable','on');
    
    % Subproject edge
    subProjectsEdge=findProjSubDir(projDir,'edge');
    set(handles.popupEdge,'String',subProjectsEdge);
    % Set as default the one chosen in the project
    [found,pos]=ismember(settings.subProjects(4),subProjectsEdge);
    if found==0
        error('The directory stored in fsmCenter''s settings does not actually exist on disk!');
    else
        set(handles.popupEdge,'Value',pos);
    end
    set(handles.popupEdge,'Enable','on');

    % Subproject lpla
    subProjectsLpla=findProjSubDir(projDir,'lpla');
    set(handles.popupLpla,'String',subProjectsLpla);
    % Set as default the one chosen in the project
    [found,pos]=ismember(settings.subProjects(2),subProjectsLpla);
    if found==0
        error('The directory stored in fsmCenter''s settings does not actually exist on disk!');
    else
        set(handles.popupLpla,'Value',pos);
    end
    set(handles.popupLpla,'Enable','on');

end

% --- Executes during object creation, after setting all properties.
function popupLpla_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupLpla (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupLpla.
function popupLpla_Callback(hObject, eventdata, handles)
% hObject    handle to popupLpla (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupLpla contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupLpla


