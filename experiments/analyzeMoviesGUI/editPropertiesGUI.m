function varargout = editPropertiesGUI(varargin)
% EDITPROPERTIESGUI M-file for editPropertiesGUI.fig
%      EDITPROPERTIESGUI, by itself, creates a new EDITPROPERTIESGUI or raises the existing
%      singleton*.
%
%      H = EDITPROPERTIESGUI returns the handle to a new EDITPROPERTIESGUI or the handle to
%      the existing singleton*.
%
%      EDITPROPERTIESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITPROPERTIESGUI.M with the given input arguments.
%
%      EDITPROPERTIESGUI('Property','Value',...) creates a new EDITPROPERTIESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before editPropertiesGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to editPropertiesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% 
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help editPropertiesGUI

% Last Modified by GUIDE v2.5 02-Jul-2003 08:43:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @editPropertiesGUI_OpeningFcn, ...
    'gui_OutputFcn',  @editPropertiesGUI_OutputFcn, ...
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


%---------------------------------------------------------------------------------------------------

% --- Executes on button press in edit_accept_all_PB.
function edit_accept_all_PB_Callback(hObject, eventdata, handles)

handles = guidata(hObject);

edit_readguidata(handles);
amgHandles = guidata(handles.amgHandles.AMG);

%get AMG to accept all
amgHandles.acceptAll = 1;
guidata(amgHandles.AMG,amgHandles);

figure(handles.amgHandles.AMG);

close(handles.EPGUI);

% --- Executes on button press in edit_apply2all_PB.
function edit_apply2all_PB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_apply2all_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

edit_readguidata(handles);
amgHandles = guidata(handles.amgHandles.AMG);

%ask if movie properties should be changed, too
ans = questdlg({'Do you want to change movie properties, too?';'(job options are changed automatically)'},'','Yes','No','Cancel','Yes');
switch ans
    case 'Yes'
        amgHandles.apply2all = 2;
    case 'No'
        amgHandles.apply2all = 1;
    otherwise
        %if cancel or close window: don't do anything
        return
end

guidata(amgHandles.AMG,amgHandles);

figure(handles.amgHandles.AMG);

close(handles.EPGUI);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%--------------------------------------------------------------------------
%-------------------Callbacks calling a function --------------------------
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = editPropertiesGUI_OutputFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_pixel_xy_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pixel_xy_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<= 0|isnan(num)
        set(hObject,'String','');
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_pixel_z_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pixel_z_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<= 0|isnan(num)
        set(hObject,'String','');
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_lensID_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_lensID_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<= 0|isnan(num)
        set(hObject,'String','');
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_wavelength_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_wavelength_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<= 0|isnan(num)
        set(hObject,'String','');
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_NA_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_NA_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<= 0|isnan(num)
        set(hObject,'String','');
    else %calculate filter parameters
        wvl = str2double(get(handles.edit_wavelength_txt,'String'));
        NA = num;
        FT_XY =  handles.sigmaCorrection(1)*(0.21*wvl/NA)/handles.pixelsizeXYZ(1);
        FT_Z =  handles.sigmaCorrection(2)*(0.66*wvl*1.33/(NA)^2)/handles.pixelsizeXYZ(2);
        patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
        handles.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
        set(handles.edit_filterX_txt,'String',sprintf('%0.3f',FT_XY));
        set(handles.edit_filterY_txt,'String',sprintf('%0.3f',FT_XY));
        set(handles.edit_filterZ_txt,'String',sprintf('%0.3f',FT_Z));
        set(handles.edit_patchX_txt,'String',patchXYZ(1));
        set(handles.edit_patchY_txt,'String',patchXYZ(2));
        set(handles.edit_patchZ_txt,'String',patchXYZ(3));
        guidata(hObject,handles);
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_timeLapse_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_timeLapse_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp = get(hObject,'String');
try %tries to convert input into double. If it fails or if input<0, input is deleted
    num = str2double(inp);
    if num<0|isnan(num)
        set(hObject,'String','');
    end
catch
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_type_PD_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in edit_type_PD.
function edit_type_PD_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_detect.
function edit_check_detect_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status); %add filter
    if ~any(marksNot2set==1)
        set(handles.checkH(1),'Value',1);
    end
    set(handles.detectH,'Enable','on');
else %unselect all not possible
    numDone = max(length(bsum2lvec(handles.status)),2);
    set(handles.checkH(numDone+1:end),'Value',0);
    set(handles.detectH,'Enable','off');
    set(handles.linkH,'Enable','off');
    %set(handles.trackH,'Enable','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_link.
function edit_check_link_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test, enable editing
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
    end
    set(handles.linkH,'Enable','on');
else %unselect all not possible, disable editing
    numDone = max(length(bsum2lvec(handles.status)),3);
    set(handles.checkH(numDone+1:end),'Value',0);
    set(handles.linkH,'Enable','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_labelgui.
function edit_check_labelgui_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
        set(handles.detectH,'Enable','on');
    end
    if ~any(marksNot2set==4) %need link
        set(handles.checkH(3),'Value',1);
        set(handles.linkH,'Enable','on');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_track.
function edit_check_track_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
        set(handles.detectH,'Enable','on');
    end
    if ~any(marksNot2set==4) %need link
        set(handles.checkH(3),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==8) %need labelgui
        set(handles.checkH(4),'Value',1);
        set(handles.linkH,'Enable','on');
    end
else %unselect all not possible
    numDone = max(length(bsum2lvec(handles.status)),5);
    set(handles.checkH(numDone+1:end),'Value',0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_labelgui2.
function edit_check_labelgui2_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
        set(handles.detectH,'Enable','on');
    end
    if ~any(marksNot2set==4) %need link
        set(handles.checkH(3),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==8) %need labelgui
        set(handles.checkH(4),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==16) %need tracker
        set(handles.checkH(5),'Value',1);
        %set(handles.trackH,'Enable','on');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_analysis1.
function edit_check_analysis1_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
        set(handles.detectH,'Enable','on');
    end
    if ~any(marksNot2set==4) %need link
        set(handles.checkH(3),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==8) %need labelgui
        set(handles.checkH(4),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==16) %need tracker
        set(handles.checkH(5),'Value',1);
        %set(handles.trackH,'Enable','on');
    end
    if ~any(marksNot2set==32) %need labelgui2
        set(handles.checkH(6),'Value',1);
        set(handles.linkH,'Enable','on');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_analysis2.
function edit_check_analysis2_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if selected, unselect test
if get(hObject,'Value')
    handles = guidata(hObject);
    marksNot2set = bsum2bvec(handles.status);
    if ~any(marksNot2set==1) %need filter
        set(handles.checkH(1),'Value',1);
    end
    if ~any(marksNot2set==2) %need detect
        set(handles.checkH(2),'Value',1);
        set(handles.detectH,'Enable','on');
    end
    if ~any(marksNot2set==4) %need link
        set(handles.checkH(3),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==8) %need labelgui
        set(handles.checkH(4),'Value',1);
        set(handles.linkH,'Enable','on');
    end
    if ~any(marksNot2set==16) %need tracker
        set(handles.checkH(5),'Value',1);
        %set(handles.trackH,'Enable','on');
    end
    if ~any(marksNot2set==32) %need labelgui2
        set(handles.checkH(6),'Value',1);
        set(handles.linkH,'Enable','on');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_filter.
function edit_check_filter_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if unselected, unselect everything else
if get(hObject,'Value')
    handles = guidata(hObject);
    set(handles.linkH,'Enable','off');
    set(handles.detectH,'Enable','off');
    if isempty(handles.header.correctInfo)
        set(handles.backgroundH,'Enable','on');
        set(handles.backgroundH(3:4),'Enable','off');
    end
    %set(handles.trackH,'Enable','off');
    %set(handles.testH,'Enable','on');
else
    set(handles.checkH(2:end),'Value',0); %do not try to assign value with a matrix: objects will disappear!
    set(handles.linkH,'Enable','off');
    set(handles.detectH,'Enable','off');
    set(handles.backgroundH,'Enable','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_OK_PB.
function edit_OK_PB_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles = guidata(hObject);

edit_readguidata(handles);

figure(handles.amgHandles.AMG);

close(handles.EPGUI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_cancel_PB.
function edit_cancel_PB_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%confirm user input
really = questdlg('Do you really want to quit without saving changes?','Quit?','yes','no','yes');
if strcmp(really,'no')
    return %end evaluation here
end

handles = guidata(hObject);
figure(handles.amgHandles.AMG);

close(handles.EPGUI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_CreateNew.
function edit_check_createNew_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if createNew is being checked, it can not be unchecked again

if get(hObject,'Value')==0;
    set(hObject,'Value',1);
    return
end

handles = guidata(hObject);
set(handles.detectH,'Enable','on');
set(handles.checkH,'Enable','on');

%if one wants to create a new file with new settings, you should not reuse
%something created with other settings (might become modal later)
handles.status = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_maxslope_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_maxslope_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_ftest_prob_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ftest_prob_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftest_prob = str2double(get(hObject,'String'));
if ~(ftest_prob>0&ftest_prob<1)
    set(hObject,'String','');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_IDopt_weight_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hints: get(hObject,'String') returns contents of edit_IDopt_weight_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_IDopt_weight_txt as a double

weight = str2double(get(hObject,'String'));
if ~(weight>=0&weight<=1)
    set(hObject,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit_check_split.
function edit_check_split_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hint: get(hObject,'Value') returns toggle state of edit_check_split
h = warndlg('Sorry, not implemented yet','Feature unavailable');
uiwait(h)
set(hObject,'Value',1-get(hObject,'Value')); %return to previous state


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================== Background package ==========================

function edit_check_correctBG_Callback(hObject, eventdata, handles)

%if being checked: activate the other BG - features
if get(hObject,'Value')
    set(handles.backgroundH([2,5]),'Enable','on');
    %if the first/last - RB is on, enable its input fields
    if get(handles.backgroundH(5),'Value')
        set(handles.backgroundH(6:end),'Enable','on');
        %change numTimepoints
        %get # of timepoints used for calibration
        noFramesL = str2double(get(handles.backgroundH(7),'String'));
        noFramesS = str2double(get(handles.backgroundH(6),'String'));
        delta = 0;
        if ~isnan(noFramesL)
            delta = delta+noFramesL;
        end
        if ~isnan(noFramesS)
            delta = delta + noFramesS;
        end
        if delta<handles.header.numTimepoints
            %change # of timepoints
            set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
                'ForegroundColor','r');
        else %we would be left with zero movie frames, so delete first/last
            set(handles.backgroundH(6),'String','');
            set(handles.backgroundH(7),'String','');
            set(handles.edit_movieLength_txt,...
                'ForegroundColor','r');
        end
    end
else %it's being unchecked: disable features and reset numTimePoints
    set(handles.backgroundH(2:end),'Enable','off');
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints),...
        'ForegroundColor','k');
    %future: reset also analyzeLimited#ofTimepoints
end
%store current BGState
writeCurrentBackgroundState(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_RB_standardBG_Callback(hObject, eventdata, handles)

%if clicked: unclick other, deactivate its input fields, activate own PDs
if get(hObject,'Value')
    set(handles.backgroundH(5),'Value',0);
    set(handles.backgroundH(6:end),'Enable','off');
    set(handles.backgroundH(3:4),'Enable','on');
    %reset # of timepoints
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints),...
        'ForegroundColor','k');
end
%store current BGState
writeCurrentBackgroundState(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_standardMovie1_PD_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit_standardMovie1_PD_Callback(hObject, eventdata, handles)

%check whether we add a new file whether we just continue with an existing one
pdString = get(hObject,'String');
currentSelection = get(hObject,'Value');

%make sure we have no wrong selection
if strcmp(pdString{currentSelection},'select movie #1')
    return
end

if strcmp(pdString{currentSelection},'add new...')
    
    %check if default biodata-dir exists
    oldDir = pwd;
    mainDir = cdBiodata(0);
    
    %find calibration movie
    [fileName,pathName] = uigetfile({'*.r3d;*.r3c','raw movie files'},'select project movie');
    
    %return if user has selected nothing
    if fileName==0 | ~strcmp(fileName(end-3:end-1),'.r3') %can be r3d or r3c
        %change back to correct dir
        cd(oldDir);
        return
    end
    
    %if last part of pathname differs from moviename: create project directory
    %and move movie and logfile; else read projectName
    calMovieName = fileName(1:end-4);
    
    %see if project is in mainDir, do it case-unsensitive
    if strcmpi(pathName(1:length(mainDir)),mainDir)==1
        relPathName = pathName(length(mainDir)+2:end-1); %begins without filesep, ends with none
    else
        ans = questdlg('moviefile is not in main file structure or you have no env-var BIODATA. There will be problems sharing your data',...
            'WARNING','Continue','Cancel','Cancel');
        if strcmp(ans,'Cancel')
            %change back to correct dir
            cd(oldDir);
            return %end evaluation here
        else
            relPathName = pathName(1:end-1); %ends without filesep
        end
    end
    
    %write filename into Pulldown-String
    pdString{end+1} = calMovieName;
    %check whether #1 is select. If yes, delete
    if strcmp(pdString{1},'select movie #1')
        pdString(1) = [];
    end
    newItem = length(pdString);
    set(hObject,'String',pdString,'Value',newItem);
    
    %store pathname&filename
    handles.correctBackground1{end+1,1} = pathName;
    handles.correctBackground1{end,2} = fileName;
    
    %store data
    guidata(hObject, handles);
    
    %change back to correct dir
    cd(oldDir);
else    
    %eliminate selectMovie-entry if necessary
    if strcmp(pdString{1},'select movie #1')
        pdString(1) = [];
        currentSelection = currentSelection -1;
    end
    set(hObject,'String',pdString,'Value',currentSelection);
end
%store current BGState
writeCurrentBackgroundState(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_standardMovie2_PD_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit_standardMovie2_PD_Callback(hObject, eventdata, handles)

%check whether we add a new file whether we just continue with an existing one
pdString = get(hObject,'String');
currentSelection = get(hObject,'Value');

%make sure we have no wrong selection
if strcmp(pdString{currentSelection},'select movie #2')
    return
end

if strcmp(pdString{currentSelection},'add new...')
    
    %check if default biodata-dir exists
    oldDir = pwd;
    mainDir = cdBiodata(0);
    
    %find calibration movie
    [fileName,pathName] = uigetfile({'*.r3d;*.r3c','raw movie files'},'select project movie');
    
    %return if user has selected nothing
    if fileName==0 | ~strcmp(fileName(end-3:end-1),'.r3') %can be r3d or r3c
        %change back to correct dir
        cd(oldDir);
        return
    end
    
    %if last part of pathname differs from moviename: create project directory
    %and move movie and logfile; else read projectName
    calMovieName = fileName(1:end-4);
    
    %see if project is in mainDir, do it case-unsensitive
    if strcmpi(pathName(1:length(mainDir)),mainDir)==1
        relPathName = pathName(length(mainDir)+2:end-1); %begins without filesep, ends with none
    else
        ans = questdlg('moviefile is not in main file structure or you have no env-var BIODATA. There will be problems sharing your data',...
            'WARNING','Continue','Cancel','Cancel');
        if strcmp(ans,'Cancel')
            %change back to correct dir
            cd(oldDir);
            return %end evaluation here
        else
            relPathName = pathName(1:end-1); %ends without filesep
        end
    end
    
    %write filename into Pulldown-String
    pdString{end+1} = calMovieName;
    %check whether #1 is select. If yes, delete
    if strcmp(pdString{1},'select movie #1')
        pdString(1) = [];
    end
    newItem = length(pdString);
    set(hObject,'String',pdString,'Value',newItem);
    
    %store pathname&filename
    handles.correctBackground2{end+1,1} = pathName;
    handles.correctBackground2{end,2} = fileName;
    
    %store data
    guidata(hObject, handles);
    
    %change back to correct dir
    cd(oldDir);
else
    %eliminate selectMovie-entry if necessary
    if strcmp(pdString{1},'select movie #2')
        pdString(1) = [];
        currentSelection = currentSelection - 1;
    end
    set(hObject,'String',pdString,'Value',currentSelection);
end
%store current BGState
writeCurrentBackgroundState(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_RB_firstAndLast_Callback(hObject, eventdata, handles)

%if clicked: unclick other, enable input fields, disable PD
if get(hObject,'Value')
    set(handles.backgroundH(2),'Value',0);
    set(handles.backgroundH(4:end),'Enable','on');
    set(handles.backgroundH(3:4),'Enable','off');
    pdString1 = get(handles.backgroundH(3),'String');
    pdString1 = [{'select movie #1'};pdString1];
    pdString2 = get(handles.backgroundH(4),'String');
    pdString2 = [{'select movie #2'};pdString2];
    set(handles.backgroundH(3:4),{'String'},{pdString1;pdString2},'Value',1);
    %change numTimepoints
    %get # of timepoints used for calibration
    noFramesL = str2double(get(handles.backgroundH(7),'String'));
    noFramesS = str2double(get(handles.backgroundH(6),'String'));
    delta = 0;
    if ~isnan(noFramesL)
        delta = delta+noFramesL;
    end
    if ~isnan(noFramesS)
        delta = delta + noFramesS;
    end
    if delta<handles.header.numTimepoints
        %change # of timepoints
        set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
            'ForegroundColor','r');
    else %we would be left with zero movie frames, so delete first/last
        set(handles.backgroundH(6),'String','');
        set(handles.backgroundH(7),'String','');
        set(handles.edit_movieLength_txt,...
            'ForegroundColor','r');
    end
end
guidata(hObject,handles);
%store current BGState
writeCurrentBackgroundState(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_txt_useLast_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit_txt_useLast_Callback(hObject, eventdata, handles)

%get # of frames for first and last. Check that we have strings and that there
%are frames left after cropping, then update numTimepoints left
noFramesL = str2double(get(handles.backgroundH(7),'String'));
noFramesS = str2double(get(handles.backgroundH(6),'String'));
delta = 0;
if ~isnan(noFramesL)&noFramesL>0&round(noFramesL)==noFramesL
    delta = delta+noFramesL;
else
    set(handles.backgroundH(7),'String','0');
end
if ~isnan(noFramesS)
    delta = delta + noFramesS;
end
if delta<handles.header.numTimepoints
    %change # of timepoints
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
        'ForegroundColor','r');
else %we would be left with zero movie frames, so delete last
    set(handles.backgroundH(7),'String','0');
    delta = delta-noFramesL;
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
        'ForegroundColor','r');
end

%store current BGState
writeCurrentBackgroundState(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_txt_useFirst_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit_txt_useFirst_Callback(hObject, eventdata, handles)

%get # of frames for first and last. Check that we have strings and that there
%are frames left after cropping, then update numTimepoints left
noFramesL = str2double(get(handles.backgroundH(7),'String'));
noFramesS = str2double(get(handles.backgroundH(6),'String'));
delta = 0;
if ~isnan(noFramesL)
    delta = delta+noFramesL;
end
if ~isnan(noFramesS)&noFramesS>0&round(noFramesS)==noFramesS
    delta = delta + noFramesS;
else
    set(handles.backgroundH(6),'String','0');
end
if delta<handles.header.numTimepoints
    %change # of timepoints
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
        'ForegroundColor','r');
else %we would be left with zero movie frames, so delete last
    set(handles.backgroundH(6),'String','0');
    delta = delta-noFramesS;
    set(handles.edit_movieLength_txt,'String',num2str(handles.header.numTimepoints-delta),...
        'ForegroundColor','r');
end

%store current BGState
writeCurrentBackgroundState(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeCurrentBackgroundState(hObject,handles)
correctBG = struct('correctFiles',[],'correctFrames',[],'header',[]);
if get(handles.edit_check_correctBG,'Value')
    %one of the two fields will not be empty
    
    correctBG.header = handles.header;
    switch get(handles.edit_RB_firstAndLast,'Value')
        case 1
            %get number of first and last dark frames
            noFramesL = str2double(get(handles.backgroundH(7),'String'));
            noFramesS = str2double(get(handles.backgroundH(6),'String'));
            if any(isnan(noFramesL)|isnan(noFramesS)|noFramesS+noFramesL<1)
                noFramesS = [];
                noFramesL = [];
            end
            correctBG.correctFrames = [noFramesS,noFramesL];
        case 0
            %get path and name of correctionFiles
            fileNum1 = get(handles.backgroundH(3),'Value');
            fileNum2 = get(handles.backgroundH(4),'Value');
            if fileNum1>1
                cfCell = handles.correctBackground1(fileNum1-1,:);
            else
                cfCell = {};
            end
            if fileNum2>1
                cfCell = [cfCell;handles.correctBackground2(fileNum2-1,:)];
            end
            correctBG.correctFiles = cfCell;
    end
else
    %if no correction, bgInfo completely empty
    correctBG = [];
end

%store correctBGInfo
handles.currentBGState = correctBG;
guidata(hObject,handles);

%--------------------------------------------------------------------------
%----------Callbacks without function--------------------------------------
%--------------------------------------------------------------------------



% --- Executes during object creation, after setting all properties.
function edit_exposureTime_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exposureTime_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_exposureTime_txt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exposureTime_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_exposureTime_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_exposureTime_txt as a double



% --- Executes during object creation, after setting all properties.
function edit_movieSize_z_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_z_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_movieSize_z_txt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_z_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_movieSize_z_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_movieSize_z_txt as a double


% --- Executes during object creation, after setting all properties.
function edit_movieSize_y_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_y_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_movieSize_y_txt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_y_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_movieSize_y_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_movieSize_y_txt as a double


% --- Executes during object creation, after setting all properties.
function edit_movieSize_x_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_x_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_movieSize_x_txt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_movieSize_x_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_movieSize_x_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_movieSize_x_txt as a double


% --- Executes during object creation, after setting all properties.
function edit_NDfilter_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NDfilter_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_NDfilter_txt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NDfilter_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NDfilter_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_NDfilter_txt as a double


% --- Executes on button press in edit_check_CC_G1.
function edit_check_CC_G1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_G1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_G1


% --- Executes on button press in edit_check_CC_SPhase.
function edit_check_CC_SPhase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_SPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_SPhase


% --- Executes on button press in edit_check_CC_Start.
function edit_check_CC_Start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_Start


% --- Executes on button press in edit_check_CC_Meta.
function edit_check_CC_Meta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_Meta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_Meta


% --- Executes on button press in edit_check_CC_G2.
function edit_check_CC_G2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_G2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_G2


% --- Executes on button press in edit_check_CC_Telo.
function edit_check_CC_Telo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_Telo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_Telo


% --- Executes on button press in edit_check_CC_Ana.
function edit_check_CC_Ana_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_CC_Ana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_CC_Ana


% --- Executes on button press in edit_check_strain_ndc80_1.
function edit_check_strain_ndc80_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_ndc80_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_ndc10_1.
function edit_check_strain_ndc10_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_ndc10_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_bik1.
function edit_check_strain_bik1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_bik1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_strain_bik1
if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end
% --- Executes on button press in edit_check_strain_ipl1_321.
function edit_check_strain_ipl1_321_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_ipl1_321 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_dam1_1.
function edit_check_strain_dam1_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_dam1_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_WT.
function edit_check_strain_WT_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_WT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(2:end),'Value',0);
end


% --- Executes on button press in edit_check_strain_ipl1_1.
function edit_check_strain_ipl1_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_ipl1_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_dam1_11.
function edit_check_strain_dam1_11_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_dam1_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_strain_bim1.
function edit_check_strain_bim1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_strain_bim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.strainH(1),'Value',0);
end

% --- Executes on button press in edit_check_drugs_alphaFactor.
function edit_check_drugs_alphaFactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_drugs_alphaFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_drugs_alphaFactor


% --- Executes on button press in edit_check_drugs_benomyl.
function edit_check_drugs_benomyl_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_drugs_benomyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_drugs_benomyl


% --- Executes on button press in edit_check_drugs_nocodazole.
function edit_check_drugs_nocodazole_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_drugs_nocodazole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_drugs_nocodazole


% --- Executes on button press in edit_check_drugs_hydroxyurea.
function edit_check_drugs_hydroxyurea_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_drugs_hydroxyurea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_drugs_hydroxyurea


% --- Executes during object creation, after setting all properties.
function edit_temperature_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_temperature_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in edit_temperature_PD.
function edit_temperature_PD_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temperature_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edit_temperature_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_temperature_PD


% --- Executes on button press in edit_check_analyzeSelected.
function edit_check_analyzeSelected_Callback(hObject, eventdata, handles)
% hObject    handle to edit_check_analyzeSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_check_analyzeSelected

function edit_startTime_txt_CreateFcn(hObject, eventdata, handles)

function edit_startTime_txt_Callback(hObject, eventdata, handles)

function edit_endTime_txt_CreateFcn(hObject, eventdata, handles)

function edit_endTime_txt_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_timeptsInMem_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_timeptsInMem_txt_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hints: get(hObject,'String') returns contents of edit_timeptsInMem_txt as text
%        str2double(get(hObject,'String')) returns contents of edit_timeptsInMem_txt as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_check_eraseAllPrevData_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hint: get(hObject,'Value') returns toggle state of edit_check_eraseAllPrevData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_IDopt_checkIntensity_PD_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in edit_IDopt_checkIntensity_PD.
function edit_IDopt_checkIntensity_PD_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hints: contents = get(hObject,'String') returns edit_IDopt_checkIntensity_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_IDopt_checkIntensity_PD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_IDopt_verbose_PD_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in edit_IDopt_verbose_PD.
function edit_IDopt_verbose_PD_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hints: contents = get(hObject,'String') returns edit_IDopt_verbose_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_IDopt_verbose_PD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_IDopt_weight_txt_CreateFcn(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



