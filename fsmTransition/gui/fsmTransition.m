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

% Last Modified by GUIDE v2.5 06-Sep-2004 14:46:16

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

function editExtractMedian_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editExtractFrames_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editExtractSupport_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSplineTolTrans_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSplineTolEdge_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSplineMedian_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editProtrusionMapsInterpY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editProtrusionMapsInterpX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function popupEdge_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function popupLpla_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editExtractMinDist_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function checkPanelSegments_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMPORT PRPANEL DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushPanel_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
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
savePath=[lplaProjDir,'prPanelDatFiles'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'prPanelDatFiles');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end
save([savePath,filesep,'normals.mat'],'normals');
save([savePath,filesep,'profiles.mat'],'profiles');
save([savePath,filesep,'protrusions.mat'],'protrusions');
save([savePath,filesep,'edgePixels.mat'],'edgePixels');
save([savePath,filesep,'segments.mat'],'segments');
msg=['Files saved in [',savePath,']'];
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

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

toggleKinetic=get(handles.checkExtractKin,'Value');
toggleSpeed=get(handles.checkExtractSpeed,'Value');

% Load profiles.mat
profilesFileName=[lplaProjDir,filesep,'prPanelDatFiles',filesep,'profiles.mat'];
if exist(profilesFileName,'file')==0
    uiwait(errordlg('Could not find file ''profiles.mat''. Make sure to ''Import prPanel data''.','Error','modal'));
    return
end
try
    s=load(profilesFileName);
    profiles=s.profiles;
catch
    uiwait(errordlg('Corrupted ''profiles.mat'' file. Try to re-import it (''Import prPanel data'') or re-run prPanel .','Error','modal'));
    return
end

% Other parameters
extractTransitionFrames=str2num(get(handles.editExtractFrames,'String'));
medianFiltFrames=str2num(get(handles.editExtractMedian,'String'));
dist=str2num(get(handles.editExtractSupport,'String'));
toggleDpRatio=get(handles.checkExtractDpRatio,'Value');
minDist=str2num(get(handles.editExtractMinDist,'String'));
DEBUG=get(handles.checkExtractDebug,'Value');

% Run function
[fPositionsScore,fPositionsVel,fPositionsTrans,dpRatio,scoreProfiles,speedProfiles]=fsmTransExtractLpToLaTransition(tackProjDir,toggleKinetic,toggleSpeed,toggleDpRatio,profiles,dist,extractTransitionFrames,DEBUG,medianFiltFrames,minDist);

% Save outputs to disk
savePath=[lplaProjDir,'lpToLaTransition'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'lpToLaTransition');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end

eval(['save ',savePath,filesep,'fPositionsScore.mat fPositionsScore']);
eval(['save ',savePath,filesep,'fPositionsVel.mat fPositionsVel']);
eval(['save ',savePath,filesep,'fPositionsTrans.mat fPositionsTrans']);
eval(['save ',savePath,filesep,'dpRatio.mat dpRatio']);
eval(['save ',savePath,filesep,'scoreProfiles.mat scoreProfiles']);
eval(['save ',savePath,filesep,'speedProfiles.mat speedProfiles']);
eval(['save ',savePath,filesep,'extractTransitionFrames.mat extractTransitionFrames']);
% Save settings
fid=fopen([savePath,filesep,'extractTransitionSettings.txt'],'a+');
if fid~=-1
    fprintf(fid,'>%s\n',datestr(now,0));
    fprintf(fid,'Use kinetic criterion              : %d\n',toggleKinetic);
    fprintf(fid,'Use speed criterion                : %d\n',toggleSpeed);
    fprintf(fid,'Calculate depoly/poly ratio        : %d\n',toggleDpRatio);
    fprintf(fid,'Median filter order                : %d\n',medianFiltFrames);
    fprintf(fid,'Support (pixels)                   : %d\n',dist);
    fprintf(fid,'Time averaging (frames)            : %d\n',extractTransitionFrames);
    fprintf(fid,'Min leading-edge distance (pixels) : %d\n\n',minDist);
    fprintf(fid,'\n');
    fclose(fid);
end
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','fPositionsScore',fPositionsScore);
assignin('base','fPositionsVel',fPositionsVel);
assignin('base','fPositionsTrans',fPositionsTrans);
assignin('base','dpRatio',dpRatio);
assignin('base','scoreProfiles',scoreProfiles);
assignin('base','speedProfiles',speedProfiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END EXTRACT TRANSITION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editExtractSupport_Callback(hObject, eventdata, handles)

function checkExtractdpRatio_Callback(hObject, eventdata, handles)

function checkExtractDpRatio_Callback(hObject, eventdata, handles)

function checkExtractDebug_Callback(hObject, eventdata, handles)
if get(handles.checkExtractDebug,'Value')==1
    uiwait(warndlg('This creates A LOT of figures. Don''t use it if you plan to analyze many frames with many profiles.','Warning','modal'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  FIT SPLINES
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushSpline_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

% Find which of the radio buttons is active ('Use ... transition')
indx=find([get(handles.radioSplineKin,'Value') get(handles.radioSplineSpeed,'Value') get(handles.radioSplineTrans,'Value')]);

switch indx
    case 1, fName=[lplaProjDir,filesep,'lpToLaTransition',filesep,'fPositionsScore.mat'];  
    case 2, fName=[lplaProjDir,filesep,'lpToLaTransition',filesep,'fPositionsSpeed.mat']; 
    case 3, fName=[lplaProjDir,filesep,'lpToLaTransition',filesep,'fPositionsTrans.mat']; 
    otherwise, error('One of the radio buttons in ''Calculate transition to edge distances'' should be on.');
end

if exist(fName,'file')==2
    s=load(fName);
    switch indx
        case 1, positionsTrans=s.fPositionsScore;
        case 2, positionsTrans=s.fPositionsSpeed;
        case 3, positionsTrans=s.fPositionsTrans;
        otherwise, error('The loaded fPositions*.mat file is not valid.');
    end
else
    uiwait(errordlg('Please run ''Extract transition'' in fsmTransition','Error','modal'));
    return
end

% Edge pixels
edgePixelsFileName=[lplaProjDir,filesep,'prPanelDatFiles',filesep,'edgePixels.mat'];
try
    e=load(edgePixelsFileName);
    edgePixels=e.edgePixels;
catch
    uiwait(errordlg('File ''edgePixels.mat'' not found or not valid. Please re-run ''Import prPanel data''.','Error','modal'));
    return
end

% Number of frames used for averaging in 'Extract Lp-to-La transition'
% This is needed to know which transition frame has to be matched with eah edge frame
% e.g. if the number of frames for averaging is 5, the transitions of frames 1 to 5 are
% averaged to deliver a transition which should then be compared with the edge at frame 3.
try
    s=load([lplaProjDir,filesep,'lpToLaTransition',filesep,'extractTransitionFrames.mat']);
    extractTransitionFrames=s.extractTransitionFrames;
catch
    uiwait(errordlg('Please run ''Extract transition'' in fsmTransition.','Error','modal'));
    return
end

% Other parameters
toleranceEdge=str2num(get(handles.editSplineTolEdge,'String'));
toleranceTrans=str2num(get(handles.editSplineTolTrans,'String'));
medianFilterFrames=str2num(get(handles.editSplineMedian,'String'));

% Run function
[spEdge_y,spEdge_x,spTrans_y,spTrans_x,splineCoordsEdge,splineCoordsTrans,numSplinePoints]=fsmTransFitSplines(positionsTrans,edgePixels,toleranceEdge,toleranceTrans,medianFilterFrames,extractTransitionFrames);

% Save output
savePath=[lplaProjDir,'splines'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'splines');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end
eval(['save ',savePath,filesep,'spEdge.mat spEdge_y spEdge_x']);
eval(['save ',savePath,filesep,'spTrans.mat spTrans_y spTrans_x']);
eval(['save ',savePath,filesep,'splineCoordsEdge.mat splineCoordsEdge']);
eval(['save ',savePath,filesep,'splineCoordsTrans.mat splineCoordsTrans']);
eval(['save ',savePath,filesep,'numSplinePoints.mat numSplinePoints']);
eval(['save ',savePath,filesep,'fitSplinesSettings.mat toleranceEdge toleranceTrans medianFilterFrames']);
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','spEdge_y',spEdge_y);
assignin('base','spEdge_x',spEdge_x);
assignin('base','spTrans_y',spTrans_y);
assignin('base','spTrans_x',spTrans_x);
assignin('base','splineCoordsEdge',splineCoordsEdge);
assignin('base','splineCoordsTrans',splineCoordsTrans);
assignin('base','numSplinePoints',numSplinePoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END OF FIT SPLINES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radioSplineKin_Callback(hObject, eventdata, handles)
set(handles.radioSplineKin,'Value',1);
set(handles.radioSplineSpeed,'Value',0);
set(handles.radioSplineTrans,'Value',0);

function radioSplineSpeed_Callback(hObject, eventdata, handles)
set(handles.radioSplineKin,'Value',0);
set(handles.radioSplineSpeed,'Value',1);
set(handles.radioSplineTrans,'Value',0);

function radioSplineTrans_Callback(hObject, eventdata, handles)
set(handles.radioSplineKin,'Value',0);
set(handles.radioSplineSpeed,'Value',0);
set(handles.radioSplineTrans,'Value',1);

function editSplineTolEdge_Callback(hObject, eventdata, handles)

function editSplineTolTrans_Callback(hObject, eventdata, handles)

function editSplineMedian_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CALCULATE TRANSITION TO LEADING EDGE DISTANCE
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushTransToEdgeDist_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

% Load spline through edge
try
    s=load([lplaProjDir,filesep,'splines',filesep,'spEdge.mat']);
    spEdge_y=s.spEdge_y;
    spEdge_x=s.spEdge_x;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Load spline through transition
try
    s=load([lplaProjDir,filesep,'splines',filesep,'spTrans.mat']);
    spTrans_y=s.spTrans_y;
    spTrans_x=s.spTrans_x;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Load number of spline points
try
    s=load([lplaProjDir,filesep,'splines',filesep,'numSplinePoints.mat']);
    numSplinePoints=s.numSplinePoints;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Get parameters from user interface
mech=get(handles.checkCalcMech,'Value');

% Run function
[transToEdgeDistNearest,transToEdgeDistMech]=fsmTransTransitionToEdgeDistances(spTrans_y,spTrans_x,spEdge_y,spEdge_x,numSplinePoints,mech);

% Save output
savePath=[lplaProjDir,'transToEdgeDistances'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'transToEdgeDistances');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end

eval(['save ',savePath,filesep,'transToEdgeDistNearest.mat transToEdgeDistNearest']);
eval(['save ',savePath,filesep,'transToEdgeDistMech.mat transToEdgeDistMech']);
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','transToEdgeDistNearest',transToEdgeDistNearest);
assignin('base','transToEdgeDistMech',transToEdgeDistMech);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END OF CALCULATE TRANSITION TO LEADING EDGE DISTANCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkCalcMech_Callback(hObject, eventdata, handles)

function editTEDTolMedian_Callback(hObject, eventdata, handles)
num=fix(str2num(get(handles.editTEDTolMedian,'String')));
if isempty(num)
    set(handles.editTEDTolMedian,'String','1');
else
    if num==0
        set(handles.editTEDTolMedian,'String','1');
    end        
end

function checkExtractKin_Callback(hObject, eventdata, handles)
if get(handles.checkExtractKin,'Value')==1
    set(handles.checkExtractDpRatio,'Enable','on');
else
    set(handles.checkExtractDpRatio,'Enable','off');
end    

function checkExtractSpeed_Callback(hObject, eventdata, handles)

function editExtractMedian_Callback(hObject, eventdata, handles)

function editExtractFrames_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CREATE/PLOT LP POLYGONS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushLpPolygons_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

% Load spline through edge (coordinate lists)
try
    s=load([lplaProjDir,filesep,'splines',filesep,'splineCoordsEdge.mat']);
    splineCoordsEdge=s.splineCoordsEdge;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Load spline through transition (coordinate lists)
try
    s=load([lplaProjDir,filesep,'splines',filesep,'splineCoordsTrans.mat']);
    splineCoordsTrans=s.splineCoordsTrans;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Plot polygons?
plotPolygons=get(handles.checkLpShowPolygons,'Value');

% Run function
lpPolygons=fsmTransCreateLpPolygons(splineCoordsEdge,splineCoordsTrans,plotPolygons);

% Save output
savePath=[lplaProjDir,'lpPolygons'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'lpPolygons');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end
eval(['save ',savePath,filesep,'lpPolygons.mat lpPolygons']);
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','lpPolygons',lpPolygons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END OF CREATE/PLOT LP POLYGONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkLpShowPolygons_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CALCULATE EDGE/TRANSITION PROTRUSIONS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushDt_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

% Load spline through edge
try
    s=load([lplaProjDir,filesep,'splines',filesep,'spEdge.mat']);
    spEdge_y=s.spEdge_y;
    spEdge_x=s.spEdge_x;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Check that at least two frames have been processed
if length(spEdge_y)<2
    uiwait(errordlg('Edge: at least two processed frames are required. Please re-run ''Extract transition'' and select more frames. If you selected only one frame in ''Edge Tracker'', make sure to re-run it, too.','Error','modal'));
    return
end

% Load spline through transition
try
    s=load([lplaProjDir,filesep,'splines',filesep,'spTrans.mat']);
    spTrans_y=s.spTrans_y;
    spTrans_x=s.spTrans_x;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Check that at least two frames have been processed
if length(spTrans_y)<2
    uiwait(errordlg('Transition: at least two processed frames are required. Please re-run ''Extract transition''.','Error','modal'));
    return
end

% Load number of knot points in splines
try
    s=load([lplaProjDir,filesep,'splines',filesep,'numSplinePoints.mat']);
    numSplinePoints=s.numSplinePoints;
catch
    uiwait(errordlg('Please run ''Fit splines'' in fsmTransition.','Error','modal'));
    return
end

% Create "spline points"
points=1:numSplinePoints;

% Check whether the mechanical model has to be used
mech=get(handles.checkProtrMech,'Value');

% Calculate transitions
[distancesNearestEdge,distancesMechEdge]=fsmTransProtrusions(spEdge_y,spEdge_x,points,mech);     % Edge
[distancesNearestTrans,distancesMechTrans]=fsmTransProtrusions(spTrans_y,spTrans_x,points,mech); % Transition

% Save output
savePath=[lplaProjDir,'protrusions'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'protrusions');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end
eval(['save ',savePath,filesep,'distancesNearestEdge.mat distancesNearestEdge']);
eval(['save ',savePath,filesep,'distancesMechEdge.mat distancesMechEdge']);
eval(['save ',savePath,filesep,'distancesNearestTrans.mat distancesNearestTrans']);
eval(['save ',savePath,filesep,'distancesMechTrans.mat distancesMechTrans']);
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','distancesNearestEdge',distancesNearestEdge);
assignin('base','distancesMechEdge',distancesMechEdge);
assignin('base','distancesNearestTrans',distancesNearestTrans);
assignin('base','distancesMechTrans',distancesMechTrans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END OF CALCULATE EDGE/TRANSITION PROTRUSIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkProtrMech_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CREATE EDGE/TRANSITION PROTRUSION MAPS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushProtrusionMaps_Callback(hObject, eventdata, handles)

% Get project information and paths
[projDir,tackProjDir,edgeProjDir,lplaProjDir]=getProjDir(handles);
if isempty(projDir)
    uiwait(errordlg('No project defined. Please create/load a project in fsmCenter.','Error','modal'));
    return
end    
if any([isempty(tackProjDir) isempty(edgeProjDir) isempty(lplaProjDir)])
    uiwait(errordlg('Could not find all (sub)projects. Quitting.','Error','modal'));
    return
end

% Try to load transToEdgeDistMech.mat (calculated with the mechanical model)
transToEdgeDist=[lplaProjDir,filesep,'transToEdgeDistances',filesep,'transToEdgeDistMech.mat'];
    
if exist(transToEdgeDist,'file')==0
    
    transToEdgeDist=[lplaProjDir,filesep,'transToEdgeDistances',filesep,'transToEdgeDistNearest.mat'];
    
    if exist(transToEdgeDist,'file')==0
        
        uiwait(errordlg('Please re-run ''Calculate transition to edge distances''.','Error','modal'));
        return
        
    end
    
end

% Load
try
    s=load(transToEdgeDist);
    fields=fieldnames(s);
    switch char(fields)
        case 'transToEdgeDistMech', transToEdgeDist=s.transToEdgeDistMech;
        case 'transToEdgeDistNearest', transToEdgeDist=s.transToEdgeDistNearest;
        otherwise, error('Wrong file.');
    end
    
catch
    uiwait(errordlg('Please re-run ''Calculate transition to edge distances''.','Error','modal'));
    return
end

% Try to load distancesMechEdge.mat (calculated with the mechanical model)
distancesEdge=[lplaProjDir,filesep,'protrusions',filesep,'distancesMechEdge.mat'];
    
if exist(distancesEdge,'file')==0
    
    distancesEdge=[lplaProjDir,filesep,'protrusions',filesep,'distancesNearestEdge.mat'];
    
    if exist(distancesEdge,'file')==0
        
        uiwait(errordlg('Please re-run ''Calculate protrusions''.','Error','modal'));
        return
        
    end
    
end

% Load
try
    s=load(distancesEdge);
    fields=fieldnames(s);
    switch char(fields)
        case 'distancesMechEdge', distancesEdge=s.distancesMechEdge;
        case 'distancesNearestEdge', distancesEdge=s.distancesNearestEdge;
        otherwise, error('Wrong file.');
    end
    
catch
    uiwait(errordlg('Please re-run ''Calculate protrusions''.','Error','modal'));
    return
end

% Try to load distancesMechTrans.mat (calculated with the mechanical model)
distancesTrans=[lplaProjDir,filesep,'protrusions',filesep,'distancesMechTrans.mat'];
    
if exist(distancesTrans,'file')==0
    
    distancesTrans=[lplaProjDir,filesep,'protrusions',filesep,'distancesNearestTrans.mat'];
    
    if exist(distancesTrans,'file')==0
        
        uiwait(errordlg('Please re-run ''Calculate protrusions''.','Error','modal'));
        return
        
    end
    
end

% Load
try
    s=load(distancesTrans);
    fields=fieldnames(s);
    switch char(fields)
        case 'distancesMechTrans', distancesTrans=s.distancesMechTrans;
        case 'distancesNearestTrans', distancesTrans=s.distancesNearestTrans;
        otherwise, error('Wrong file.');
    end
    
catch
    uiwait(errordlg('Please re-run ''Calculate protrusions''.','Error','modal'));
    return
end

% Other parameters
interpStepY=str2num(get(handles.editProtrusionMapsInterpY,'String'));
interpStepX=str2num(get(handles.editProtrusionMapsInterpX,'String'));
[transProtrusionMap2D,edgeProtrusionMap2D]=fsmTransProtrusionMaps(transToEdgeDist,distancesTrans,distancesEdge,[interpStepY interpStepX]);

% Save output
savePath=[lplaProjDir,'protrusions'];
if isdir(savePath)==0
    success=mkdir(lplaProjDir,'protrusions');
    if success==0
        errmsg=['Could not write to ',lplaProjDir];
        uiwait(errordlg(errmsg,'Error','modal'));
        return
    end
end
eval(['save ',savePath,filesep,'transProtrusionMap2D.mat transProtrusionMap2D']);
eval(['save ',savePath,filesep,'edgeProtrusionMap2D.mat edgeProtrusionMap2D']);
msg=['Files saved in [',savePath,']'];
uiwait(msgbox(msg,'Info','modal'));

% Return to Matlab workspace
assignin('base','transProtrusionMap2D',transProtrusionMap2D);
assignin('base','edgeProtrusionMap2D',edgeProtrusionMap2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END OF CREATE EDGE/TRANSITION PROTRUSION MAPS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editProtrusionMapsInterpX_Callback(hObject, eventdata, handles)

function editProtrusionMapsInterpY_Callback(hObject, eventdata, handles)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ACCESSORY CALLBACKS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function popupEdge_Callback(hObject, eventdata, handles)

function popupTack_Callback(hObject, eventdata, handles)

function popupLpla_Callback(hObject, eventdata, handles)

% Menu callback
function menuExit_Callback(hObject, eventdata, handles)

% Close request function
function fsmTransition_CloseRequestFcn(hObject, eventdata, handles)
fsmH=findall(0,'Tag','fsmTransition','Name','fsmTransition'); % Get the handle of fsmPostProc

% Check whether fsmPostProc is calling its own closeRequestFcn or whether fsmCenter is doing it
hFsmCenter=findall(0,'Tag','fsmCenter','Name','fsmCenter');

if hFsmCenter==hObject
    % Yes, it is fsmCenter
    force=1;
else
    % No, someone else is trying to close it
    force=0;
end

% Close
if fsmH~=hObject & force==1
    delete(fsmH); % Force closing
else
    choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
    switch choice,
        case 'Yes', delete(fsmH);
        case 'No', return;
    end
end


