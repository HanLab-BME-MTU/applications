function varargout = make_adhesion_mask(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @make_adhesion_mask_OpeningFcn, ...
                   'gui_OutputFcn',  @make_adhesion_mask_OutputFcn, ...
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


% --- Executes just before make_adhesion_mask is made visible.
function make_adhesion_mask_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to make_adhesion_mask (see VARARGIN)

% Choose default command line output for make_adhesion_mask
handles.output = hObject;
handles.data.imfile_name = [];
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = make_adhesion_mask_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cutoff_text_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_text as text
%        str2double(get(hObject,'String')) returns contents of cutoff_text as a double


% --- Executes during object creation, after setting all properties.
function cutoff_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function sigma_text_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_text as text
%        str2double(get(hObject,'String')) returns contents of sigma_text as a double


% --- Executes during object creation, after setting all properties.
function sigma_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function winsize_text_Callback(hObject, eventdata, handles)
% hObject    handle to winsize_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winsize_text as text
%        str2double(get(hObject,'String')) returns contents of winsize_text as a double


% --- Executes during object creation, after setting all properties.
function winsize_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winsize_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function thresh_text_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh_text as text
%        str2double(get(hObject,'String')) returns contents of thresh_text as a double


% --- Executes during object creation, after setting all properties.
function thresh_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function imfile_text_Callback(hObject, eventdata, handles)
% hObject    handle to imfile_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imfile_text as text
%        str2double(get(hObject,'String')) returns contents of imfile_text as a double


% --- Executes during object creation, after setting all properties.
function imfile_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imfile_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in imfile_browse.
function imfile_browse_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.tif', 'Select a .tif image file.');
if ~isequal(filename,0) 
    handles.data.imfile_name = fullfile(pathname, filename);
    if size(handles.data.imfile_name,2) > 55
        disp_string = ['... ',handles.data.imfile_name(end-55:end)];
    else
        disp_string = handles.data.imfile_name;
    end
    set(handles.imfile_text,'string',disp_string);
end
guidata(hObject, handles);


% --- Executes on button press in preview_button.
function preview_button_Callback(hObject, eventdata, handles)
if isempty(handles.data.imfile_name)
    errordlg('Please specify an image file.','Error');
    return;
end
    
written = str2double(get(handles.cutoff_text,'String'));
if ~isnan(written) && written >= 0
    cutoff = round(written);
else
    errordlg('The cutoff must be a positive integer or zero.','Error');
    return;
end

written = str2double(get(handles.sigma_text,'String'));
if ~isnan(written) && written >= 0
    sigma = written;
else
    errordlg('Sigma must be positive or zero.','Error');
    return;
end

written = str2double(get(handles.winsize_text,'String'));
if ~isnan(written) && written >= 0
    filter_mask_size = round(written);
else
    errordlg('The windowsize must be a positive integer or zero.','Error');
    return;
end

written = str2double(get(handles.thresh_text,'String'));
if ~isnan(written) && written >= 0 && written <= 1
    thresh = written;
else
    errordlg('The brightness threshold must be in the interval [0,1].','Error');
    return;
end

bild0 = imread(handles.data.imfile_name);
bild0 = (bild0 -min(min(bild0))+1);
filt_win = fspecial('gaussian',cutoff,sigma);
  
mean_kernel = double(1/(filter_mask_size*filter_mask_size)*ones(filter_mask_size,filter_mask_size));
mean_bild0 = conv2(double(ceil(bild0)), mean_kernel,'same');
bild = uint16(ceil((double(bild0) - mean_bild0 - min(min((double(bild0) - mean_bild0))) +1)));
bbw = filter2(filt_win,bild);
adh_mask = (bbw/(max(bbw(:))) >= thresh);


adh_mask(1:filter_mask_size +3,:) = 0;
adh_mask(:, 1:filter_mask_size +3) = 0;
adh_mask(:,end-filter_mask_size-3:end) = 0;
adh_mask(end-filter_mask_size-3:end,:) = 0;

axes(handles.axes1);
cla; axis equal, hold on; colormap gray, imagesc(bild0); 
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[])
hold off;

axes(handles.axes2);
cla; axis equal, hold on; colormap gray, imagesc(adh_mask);
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[])
hold off;

handles.data.adh_mask = adh_mask;
set(handles.save_button,'Enable','on');
guidata(hObject, handles);

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
if ~isempty(handles.data.adh_mask)
    adh_mask = handles.data.adh_mask;
    [vert_pos, hor_pos] = meshgrid(1:size(adh_mask,2), 1:size(adh_mask,1));
    save([handles.data.imfile_name(1:end-4),'_Adh_Mask.mat'],'adh_mask','vert_pos','hor_pos');
    set(handles.save_button,'Enable','off');
    set(handles.imfile_text,'string','');
    guidata(hObject, handles);
    axes(handles.axes1);cla;
    axes(handles.axes2);cla;
else
    errordlg('Please first take a look at the result. Not saving!','Error');
    return;
end

% --- Executes on button press in quit_button.
function quit_button_Callback(hObject, eventdata, handles)
delete(handles.make_adhesion_mask);
return;



