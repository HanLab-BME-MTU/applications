% This program (regularized fourier transform traction force
% reconstruction) was produced at the University of Heidelberg, BIOMS
% group of Ulich Schwarz. It calculates traction from a gel displacement
% field.
%
% Benedikt Sabass 20-5-2007

function varargout = preview_settings(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preview_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @preview_settings_OutputFcn, ...
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


% --- Executes just before preview_settings is made visible.
function preview_settings_OpeningFcn(hObject, eventdata, handles, varargin)

handles.data = varargin{1};

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
set(handles.preview_image,'String',sorted_names,'Value',1)
for i = 1:length(handles.data.strain)
    hilf_cell{i,1} = num2str(i);
end
set(handles.preview_frame,'String',hilf_cell,'Value',1);


max_eck(1:2) = [max(handles.data.strain(1).pos(:,1)), max(handles.data.strain(1).pos(:,2))];
min_eck(1:2) = [min(handles.data.strain(1).pos(:,1)), min(handles.data.strain(1).pos(:,2))];
  
handles.data.meshsize = round(sqrt((max_eck(1)-min_eck(1))*(max_eck(2)-min_eck(2))/size(handles.data.strain(1).pos,1)));
set(handles.meshsize_win,'String',num2str(handles.data.meshsize));

handles.data.regparam = 0.2/handles.data.young;
set(handles.regparam_win,'String',num2str(handles.data.regparam,'%f'));

handles.data.pos = 0;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = preview_settings_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



function meshsize_win_Callback(hObject, eventdata, handles)

written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.meshsize = written;
else
    errordlg('The mesh size must be given as a positive number.','Error');
    set(handles.meshsize_win,'String', num2str(handles.data.meshsize));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function meshsize_win_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function regparam_win_Callback(hObject, eventdata, handles)

written = str2double(get(hObject,'String'));
if ~isnan(written) && written >= 0
    handles.data.regparam = written;
else
    errordlg('The smoothing parameter must be larger or equal zero.','Error');
    set(handles.regparam_win,'String', num2str(handles.data.regparam));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function regparam_win_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in preview_image.
function preview_image_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function preview_image_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in preview_frame.
function preview_frame_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function preview_frame_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)

frame = get(handles.preview_frame,'value');
if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
    errordlg('The selected dataset is empty or not valid.','Error');
    return;
end

bild_datei_index = get(handles.preview_image,'value');
bild_dateien = get(handles.preview_image,'string');
bild_datei = bild_dateien{bild_datei_index};

[grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.strain(frame).pos, handles.data.strain(frame).vec,handles.data.meshsize);

%achim added:
figure(frame)
quiver(handles.data.strain(frame).pos(:,1),handles.data.strain(frame).pos(:,2),handles.data.strain(frame).vec(:,1),handles.data.strain(frame).vec(:,2))
title('displacement field')

%achim added: young modulus is entered in Pa NOT kPa!!!
[handles.data.pos,handles.data.vec,handles.data.force,schrott1,schrott2,handles.data.f_mat] = reg_fourier_TFM(grid_mat,u,handles.data.young,handles.data.poisson, handles.data.pix_durch_my, handles.data.meshsize, i_max,j_max,handles.data.regparam);
handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));

axes(handles.axes1);
cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild);

hilf = get(handles.show_vectors,'Value');
if hilf
    quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
end
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[])
hold off;

axes(handles.axes2);
cla; hold on; colormap jet;
fnorm = (handles.data.f_mat(:,:,2).^2 + handles.data.f_mat(:,:,1).^2).^0.5;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[]),hold off;
guidata(hObject, handles);

% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)

delete(handles.traction_reconstruction_figure);
return;


% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)

grid_mat = [];

haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
framenumber = length(handles.data.strain);
for frame = 1:framenumber
    if isempty(handles.data.strain(frame).pos) || isempty(handles.data.strain(frame).vec)|| length(handles.data.strain(frame).pos) < 2
        disp(['There dataset of frame ', num2str(frame),' was not valid. Skipping this frame!']);
        continue;
    end
    waitbar(frame/framenumber);
    [grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.strain(frame).pos, handles.data.strain(frame).vec,handles.data.meshsize, grid_mat);
    [TFM_results(frame).pos,TFM_results(frame).vec,TFM_results(frame).traction,TFM_results(frame).traction_magnitude,TFM_results(frame).energy] = reg_fourier_TFM(grid_mat,u,handles.data.young,handles.data.poisson, handles.data.pix_durch_my, handles.data.meshsize, i_max,j_max,handles.data.regparam);
end

TFM_settings.poisson = handles.data.poisson;
TFM_settings.young = handles.data.young;
TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
TFM_settings.regularization_param = handles.data.regparam;
TFM_settings.meshsize = handles.data.meshsize;
close(haha);

savefile_name = fullfile(handles.data.targetdir_name,['Reg-FTTC_results_',datestr(now, 'dd-mm-yy'),'.mat']);
if exist(savefile_name)
    button = questdlg('The file exists already. Overwrite?','Error','Yes');
    if strcmpi(button,'No') || strcmpi(button,'')
            return;
    end
end
save(savefile_name,'TFM_results','TFM_settings','-mat');

if ~isempty(handles.data.maskdir_name)
    for frame = 1:framenumber
        if exist(fullfile(handles.data.maskdir_name,['mask_',num2str(frame),'.jpg']),'file') && ~isempty(TFM_results(frame).pos) && ~isempty(TFM_results(frame).traction)
            disp(['Evaluating mask for frame: ', num2str(frame)]);
            mask = imread(fullfile(handles.data.maskdir_name,['mask_',num2str(frame),'.jpg']));
            mask = im2bw(mask);
            [vert_pos, hor_pos] = meshgrid(1:size(mask,2), 1:size(mask,1));
            C = contourc(double(mask),1);
            
            n = 1;
            i = 1;
            sites = [];
            while i < size(C,2)
                sites(n).contour(1:C(2,i),1:2) = C(1:2,i+(1:C(2,i)))';
                sites(n).contour(1:C(2,i),1) = sites(n).contour(1:C(2,i),1);
                sites(n).contour(1:C(2,i),2) = sites(n).contour(1:C(2,i),2);
                i = i+C(2,i)+1;
                n = n +1;
            end
            
            hhh = figure('Visible','on'); hold on, 
            for i=1:size(sites,2)
                plot(sites(i).contour(:,1),sites(i).contour(:,2),'b');
                text(mean(sites(i).contour(:,1)),mean(sites(i).contour(:,2)), num2str(i),'Color','k');
                inpos = inpolygon(TFM_results(frame).pos(:,1)+ TFM_results(frame).vec(:,1),TFM_results(frame).pos(:,2) + TFM_results(frame).vec(:,2),sites(i).contour(:,1),sites(i).contour(:,2));
                traction_data(i).pos(1:2) =  [mean(sites(i).contour(:,1)),mean(sites(i).contour(:,2))];
                
                if nnz(inpos) > 0
                    traction_data(i).maxabs = max((TFM_results(frame).traction(inpos,1).^2 + TFM_results(frame).traction(inpos,2).^2).^0.5);
                    traction_data(i).meantraction(1:2) = [mean(TFM_results(frame).traction(inpos,1)),mean(TFM_results(frame).traction(inpos,2))];
                    traction_data(i).meanabs = norm(traction_data(i).meantraction(1:2));
                  
                else
                    traction_data(i).maxabs = nan;
                    traction_data(i).meantraction(1:2) = nan;
                    traction_data(i).meanabs = nan;
                end
    
                sites(i).area = polyarea(sites(i).contour(:,1),sites(i).contour(:,2));
            end
            hold off;
            saveas(hhh,fullfile(handles.data.targetdir_name,['mask_locations_',num2str(frame),'.fig']),'fig')
            close(hhh);
            save(fullfile(handles.data.targetdir_name,['mask_data_',num2str(frame),'.mat']), 'sites','traction_data','-mat');
            
        end
    end
end
msgbox(['Traction data is now saved in ',' Reg-FTTC_results_',datestr(now, 'dd-mm-yy'),'.mat '],'Allesinbutter');
delete(handles.traction_reconstruction_figure);
return;


% --- Executes on button press in show_vectors.
function show_vectors_Callback(hObject, eventdata, handles)

if handles.data.pos ~=0
    axes(handles.axes1);
    cla; axis equal, hold on;  imagesc(handles.data.bild);
    if (get(hObject,'Value') == get(hObject,'Max'))
        quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
    end
    hold off;
end


