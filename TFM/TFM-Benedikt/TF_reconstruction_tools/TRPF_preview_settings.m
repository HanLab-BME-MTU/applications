function varargout = TRPF_preview_settings(varargin)
% TRPF_PREVIEW_SETTINGS M-file for TRPF_preview_settings.fig
%      TRPF_PREVIEW_SETTINGS, by itself, creates a new TRPF_PREVIEW_SETTINGS or raises the existing
%      singleton*.
%
%      H = TRPF_PREVIEW_SETTINGS returns the handle to a new TRPF_PREVIEW_SETTINGS or the handle to
%      the existing singleton*.
%
%      TRPF_PREVIEW_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRPF_PREVIEW_SETTINGS.M with the given input arguments.
%
%      TRPF_PREVIEW_SETTINGS('Property','Value',...) creates a new TRPF_PREVIEW_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TRPF_preview_settings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TRPF_preview_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TRPF_preview_settings

% Last Modified by GUIDE v2.5 17-Jun-2007 12:20:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TRPF_preview_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @TRPF_preview_settings_OutputFcn, ...
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


% --- Executes just before TRPF_preview_settings is made visible.
function TRPF_preview_settings_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;


handles.data = varargin{1};

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
set(handles.image_popup,'String',sorted_names,'Value',1)
for i = 1:length(handles.data.strain)
    hilf_cell{i,1} = num2str(i);
end
set(handles.frame_popup,'String',hilf_cell,'Value',1);

if ~isempty(sorted_names)
    bild_datei = sorted_names{1};
    handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));

    axes(handles.axes1);
    cla; axis equal, colormap gray, imagesc(handles.data.bild); 
end


max_eck(1:2) = [max(handles.data.strain(1).pos(:,1)), max(handles.data.strain(1).pos(:,2))];
min_eck(1:2) = [min(handles.data.strain(1).pos(:,1)), min(handles.data.strain(1).pos(:,2))];
  
handles.data.cutoff = 1.5*round(sqrt((max_eck(1)-min_eck(1))*(max_eck(2)-min_eck(2))/size(handles.data.strain(1).pos,1)));
set(handles.cutoff_text,'String',num2str(handles.data.cutoff));

handles.data.regparam = 0.03;
set(handles.regparam_text,'String',num2str(handles.data.regparam,'%f'));

handles.data.forcepos = [];
handles.data.redo_data = true;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TRPF_preview_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TRPF_preview_settings_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on selection change in image_popup.
function image_popup_Callback(hObject, eventdata, handles)

bild_datei_index = get(handles.image_popup,'value');
bild_dateien = get(handles.image_popup,'string');
bild_datei = bild_dateien{bild_datei_index};
handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));

axes(handles.axes1);
cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild); 
if ~isempty(handles.data.forcepos) 
    scatter(handles.data.forcepos(:,1),handles.data.forcepos(:,2),'yo');
end
hold off;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function image_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in frame_popup.
function frame_popup_Callback(hObject, eventdata, handles)
handles.data.redo_data = true;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function frame_popup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preview_button.
function preview_button_Callback(hObject, eventdata, handles)
    if isempty(handles.data.forcepos)
        errordlg('Please select locations for the Point Forces in the left image first.','Error');
        return;
    end
    frame = get(handles.frame_popup,'value');
    if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
        errordlg('The selected dataset is empty or not valid.','Error');
        return;
    end
    
    if handles.data.redo_data
        hh = waitbar(1,'Please wait','WindowStyle','modal');
        [handles.data.G,handles.data.U,handles.data.s,handles.data.V,handles.data.u] = make_data(frame, handles);
        close(hh);
        handles.data.redo_data = false;
    end
    [F,rho,eta]   = tikhonov(handles.data.U,handles.data.s,handles.data.V,handles.data.u,handles.data.regparam);
    
    forces(1:length(F)/2,1) = F(1:2:end); 
    forces(1:length(F)/2,2) = F(2:2:end); 
    
    axes(handles.axes1);
    cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild); hold on;
    quiver(handles.data.strain(frame).pos(:,1), handles.data.strain(frame).pos(:,2), handles.data.strain(frame).vec(:,1), handles.data.strain(frame).vec(:,2),'g');
    quiver(handles.data.forcepos(:,1),handles.data.forcepos(:,2), forces(:,1), forces(:,2),'r');
    hold off;
    
    guidata(hObject, handles);

% --- Executes on button press in lcurve_button.
function lcurve_button_Callback(hObject, eventdata, handles)
    if isempty(handles.data.forcepos)
        errordlg('Please select locations for the Point Forces in the left image first.','Error');
        return;
    end
    frame = get(handles.frame_popup,'value');
    if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
        errordlg('The selected dataset is empty or not valid.','Error');
        return;
    end
    axes(handles.axes2);
    cla, title('Busy');
    if handles.data.redo_data
        hh = waitbar(1,'Please wait','WindowStyle','modal');
        [handles.data.G,handles.data.U,handles.data.s,handles.data.V,handles.data.u] = make_data(frame, handles);
        close(hh);
        handles.data.redo_data = false;
    end
    [reg_corner,rho,eta,reg_param] = l_curve(handles.data.U,handles.data.s,handles.data.u,'Tikh');
    set(gca, 'YLim',[min(eta),max(eta)]),hold off;
    
    guidata(hObject, handles);


function regparam_text_Callback(hObject, eventdata, handles)

written = str2double(get(hObject,'String'));
if ~isnan(written) && written >= 0
    handles.data.regparam = written;
else
    errordlg('The regulariztion parameter must be larger or equal zero.','Error');
    set(handles.regparam_text,'String', num2str(handles.data.regparam));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function regparam_text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cutoff_text_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.cutoff = written;
else
    errordlg('The cutoff must be larger than zero.','Error');
    set(handles.cutoff_text,'String', num2str(handles.data.cutoff));
end
guidata(hObject, handles);

handles.data.redo_data = true;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cutoff_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in sframe_button.
function sframe_button_Callback(hObject, eventdata, handles)
   if isempty(handles.data.forcepos)
        errordlg('Please select locations for the Point Forces in the left image first.','Error');
        return;
    end
    frame = get(handles.frame_popup,'value');
    if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
        errordlg('The selected dataset is empty or not valid.','Error');
        return;
    end
    
    if handles.data.redo_data
        hh = waitbar(1,'Please wait','WindowStyle','modal');
        [handles.data.G,handles.data.U,handles.data.s,handles.data.V,handles.data.u] = make_data(frame, handles);
        close(hh);
        handles.data.redo_data = false;
    end
    [F,rho,eta]   = tikhonov(handles.data.U,handles.data.s,handles.data.V,handles.data.u,handles.data.regparam);
    
    forces(1:length(F)/2,1) = F(1:2:end); 
    forces(1:length(F)/2,2) = F(2:2:end); 
    
    axes(handles.axes1);
    cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild); hold on;
    quiver(handles.data.strain(frame).pos(:,1), handles.data.strain(frame).pos(:,2), handles.data.strain(frame).vec(:,1), handles.data.strain(frame).vec(:,2),'g');
    quiver(handles.data.forcepos(:,1),handles.data.forcepos(:,2), forces(:,1), forces(:,2),'r');
    hold off;
    
    
    TFM_results.pos = handles.data.strain(frame).pos;
    TFM_results.vec = handles.data.strain(frame).vec;
    TFM_results.forcepos = handles.data.forcepos;
    TFM_results(frame).forces = (handles.data.pix_durch_my^2).*handles.data.young/1000.*forces; 
   

    TFM_settings.poisson = handles.data.poisson;
    TFM_settings.young = handles.data.young;
    TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
    TFM_settings.regularization_param = handles.data.regparam;
    TFM_settings.cutoff = handles.data.cutoff;


    savefile_name = fullfile(handles.data.targetdir_name,['TRPF_results_',datestr(now, 'dd-mm-yy'),'_frame_',num2str(frame),'.mat']);
    if exist(savefile_name)
        button = questdlg('The file exists already. Overwrite?','Error','Yes');
        if strcmpi(button,'No') || strcmpi(button,'')
                return;
        end
    end
    save(savefile_name,'TFM_results','TFM_settings','-mat');
    guidata(hObject, handles);


% --- Executes on button press in analyze_button.
function analyze_button_Callback(hObject, eventdata, handles)
if isempty(handles.data.forcepos)
    errordlg('Please select locations for the Point Forces in the left image first.','Error');
    return;
end
    
haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
framenumber = length(handles.data.strain);

TFM_results(:).pos = handles.data.strain(:).pos;
TFM_results(:).vec = handles.data.strain(:).vec;
for frame = 1:framenumber
    if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
        disp(['There dataset of frame ', num2str(frame),' was not valid. Skipping this frame!']);
        continue;
    end
    waitbar(frame/framenumber);
  
    [G,U,s,V,u] = make_data(frame, handles);  
    [F,rho,eta]   = tikhonov(U,s,V,u,handles.data.regparam);
    TFM_results(frame).forcepos = handles.data.forcepos;
    TFM_results(frame).forces(1:length(F)/2,1) = (handles.data.pix_durch_my^2).*handles.data.young/1000.*F(1:2:end); 
    TFM_results(frame).forces(1:length(F)/2,2) = (handles.data.pix_durch_my^2).*handles.data.young/1000.*F(2:2:end); 
end

TFM_settings.poisson = handles.data.poisson;
TFM_settings.young = handles.data.young;
TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
TFM_settings.regularization_param = handles.data.regparam;
TFM_settings.cutoff = handles.data.cutoff;
close(haha);

savefile_name = fullfile(handles.data.targetdir_name,['TRPF_results_',datestr(now, 'dd-mm-yy'),'.mat']);
if exist(savefile_name)
    button = questdlg('The file exists already. Overwrite?','Error','Yes');
    if strcmpi(button,'No') || strcmpi(button,'')
            return;
    end
end
save(savefile_name,'TFM_results','TFM_settings','-mat');
msgbox(['Data is now saved in ','- TRPF_results_',datestr(now, 'dd-mm-yy'),'.mat -'],'Allesinbutter');
delete(handles.TRPF_preview_settings);
return;


function ppf_button_Callback(hObject, eventdata, handles)
axes(handles.axes1);
axis equal;  colormap gray; imagesc(handles.data.bild);hold on; 
if ~isempty(handles.data.forcepos)
    scatter(handles.data.forcepos(:,1),handles.data.forcepos(:,2),'yo');
end
handles.data.redo_data = true;
while 1
    [x,y,b] = ginput(1);
    if isempty(x)
        break;
    elseif b == 1
        scatter(x,y,'yo');
        handles.data.forcepos(end+1,1:2) = [x,y]';
    else
        dist = ((handles.data.forcepos(:,1) - x).^2 +(handles.data.forcepos(:,2) - y).^2).^0.5;
        [schrott,i] = min(dist);
        handles.data.forcepos(i,:) = [];
        hold off;  axis equal;  colormap gray; imagesc(handles.data.bild); hold on;
        if ~isempty(handles.data.forcepos)
            scatter(handles.data.forcepos(:,1),handles.data.forcepos(:,2),'yo');
        end
    end
end
hold off;
guidata(hObject, handles);


function [G,U,s,V,u] = make_data(frame, handles)
   
    for i = 1:size(handles.data.strain(frame).pos,1)
       for j = 1:size(handles.data.forcepos,1)
          G(2*i-1:2*i,2*j-1:2*j) = greens(handles.data.strain(frame).pos(i,1:2)-handles.data.forcepos(j,1:2),handles.data.cutoff);
       end
    end
    [U,s,V] = csvd(G);

    u(1:2:size(handles.data.strain(frame).vec,1)*2,1) = handles.data.strain(frame).vec(:,1);
    u(2:2:size(handles.data.strain(frame).vec,1)*2,1) = handles.data.strain(frame).vec(:,2);
   

function G_component = greens(r,cutoff)

    x = r(1);
    y = r(2);
    mag = sqrt(x^2+y^2);

    if (mag <= cutoff)
        G_component(1:2,1:2) = 0;
    else

        G_component(1,1) = 3/(4*pi*mag)*(1+x*x/mag/mag);
        G_component(1,2) = 3*x*y/(4*pi*mag*mag*mag);
        G_component(2,1) = G_component(1,2);
        G_component(2,2) = 3/(4*pi*mag)*(1+y*y/mag/mag);
    end


    
