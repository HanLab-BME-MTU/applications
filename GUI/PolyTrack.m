function varargout = PolyTrack(varargin)
% GUI_START M-file for GUI_start.fig
%      GUI_START, by itself, creates a new GUI_START or raises the existing
%      singleton*.
%
%      H = GUI_START returns the handle to a new GUI_START or the handle to
%      the existing singleton*.
%
%      GUI_START('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_START.M with the given input arguments.
%
%      GUI_START('Property','Value',...) creates a new GUI_START or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_start_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_start_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_start

% Last Modified by GUIDE v2.5 23-Feb-2004 14:33:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_start_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_start_OutputFcn, ...
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








% --- Executes just before GUI_start is made visible.
function GUI_start_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_start (see VARARGIN)

% Choose default command line output for GUI_start
handles.output = hObject;

defaultjob = struct('imagedirectory',[],'imagename',[],'firstimage',1,'lastimage',[],...
                   'increment',1,'savedirectory',[],'fi_background',[],'fi_nucleus',[],...
                   'la_background',[],'la_nucleus',[],'maxsearch',48,'mmpixel',[],...
                   'minsize',300,'maxsize',1500,'minsdist',30,'fi_halolevel',[],'la_halolevel',[],...
                   'minedge',10,'sizetemplate',41,'boxsize',141,'noiseparameter',0.15,...
                   'mincorrqualtempl',0.2,'leveladjust',0.7,'timestepslide',5,'mintrackcorrqual',0.5,...
                   'coordinatespicone',[],'intensityMax',4095,'bitdepth',12,'bodyname',[],'imagenameslist',[],'timeperframe',[]) ;
               
handles.defaultjob = defaultjob;

set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_start wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = GUI_start_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_job_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_job_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_job_lb.
function GUI_st_job_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_job_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_job_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_job_lb
handles = guidata(hObject);


projNum = get(hObject,'Value');

fillFields(handles,handles.jobs(projNum))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in GUI_st_addjob_pb.
function GUI_st_addjob_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_addjob_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

listhandle = handles.GUI_st_job_lb;

[filename,imagedirectory] = uigetfile({'*.tif','TIFF-files';'jobvalues.mat','Saved Values';'*.*','all files'},...
                                     'Please select an image file or jobvalues.mat');

if filename == 0
    return
end

gotvals = 0
if strcmp(filename,'jobvalues.mat')
    cd(imagedirectory)
    load jobvalues.mat
    filename = jobvalues.imagename;
    gotvals = 1;
end



%write project name into listbox
jobList = get(handles.GUI_st_job_lb,'String');
if ~iscell(jobList)
    jobList = cellstr(filename);
else
    jobList(end+1) = cellstr(filename);
end
set(handles.GUI_st_job_lb,'String',jobList);


%number of already selected projects
projNum = length(jobList); 


if gotvals==1
    handles.jobs(projNum) = jobvalues;
    clear jobvalues
else

	handles.jobs(projNum) = handles.defaultjob;
	handles.jobs(projNum).imagedirectory = imagedirectory;
	handles.jobs(projNum).imagename = filename;
	
    number = 0;
    countNum = 0;
    while ~isnan(number)
         countNum = countNum+1
         number = str2num(filename(end-(4+countNum):end-4));
         
    end
        
    
    handles.jobs(projNum).bodyname = filename(1:(end-(4+countNum)));
     
	%select current project
	set(handles.GUI_st_job_lb,'Value',projNum);
	
	
	dirList = dir(imagedirectory);
	
    dirList  =  struct2cell(dirList);
    dirList = dirList(1,:);
    ind = strmatch(handles.jobs(projNum).bodyname,dirList);
    dirList = dirList(ind)';
    handles.jobs(projNum).lastimage = length(dirList);
      
    for jRearange = 1:length(dirList)
        tmpName = char(dirList(jRearange));
        imageNum(jRearange) = str2num(tmpName(length(handles.jobs(projNum).bodyname)+1:end-4));
    end
    
    [junk,indVec] = sort(imageNum);
    
    handles.jobs(projNum).imagenameslist = dirList(indVec);
    
    
    
    
% 	howmuchindir=length(dirList);
% 	
% 	ct=1;
% 	while  ct < howmuchindir
%         if ~dirList(ct).isdir & (length(dirList(ct).name)>6)
%             if strcmp (dirList(ct).name(end-3:end),'.tif') 
%                handles.jobs(projNum).firstimage=str2num(dirList(ct).name(end-6:end-4));
%                 ct=howmuchindir+100;
%             end
%         end
%             ct=ct+1;
% 	end
% 	
% 	ct= howmuchindir;
% 	while ct >0
%         if ~dirList(ct).isdir & (length(dirList(ct).name)>6)
%             if strcmp (dirList(ct).name(end-3:end),'.tif') 
%                 handles.jobs(projNum).lastimage=str2num(dirList(ct).name(end-6:end-4));
%                 ct=-77;
%             end
%         end
%      
%         ct=ct-1;
%         
% 	end
	
    
    
    
	cd(imagedirectory)
	done = 0;
	counter = 1
	while done==0
           newdirname = [];
           newdirname = ['results',filename,num2str(counter)];
        
          if exist(newdirname,'dir')==0
             mkdir(imagedirectory,newdirname);
             tempname = [imagedirectory,newdirname];
             mkdir(tempname,'body');
             
             handles.jobs(projNum).savedirectory = [imagedirectory, newdirname];
             done = 1;
         end
         counter = counter+1;
     end
     
 end




% Update handles structure
guidata(hObject, handles);
handles = guidata(hObject);

fillFields(handles,handles.jobs(projNum))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in GUI_st_deletejob_pb.
function GUI_st_deletejob_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_deletejob_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
none = [];
handles = guidata(hObject);

jobList = get(handles.GUI_st_job_lb,'String');
projNum = get(handles.GUI_st_job_lb,'Value');
if ~iscell(jobList)
    return
    %no cell, no job, no delete
    disp('Why delete when there is nothing to delete?');
end

if length(jobList)==1
    jobList = char('No project loaded');
    none = 1
else
    jobList(projNum) = [];
end
%store new jobList
set(handles.GUI_st_job_lb,'Value',1);
set(handles.GUI_st_job_lb,'String',jobList);
%store job data
handles.jobs(projNum) = [];

guidata(hObject,handles);

if isempty(none)
    fillFields(handles,handles.jobs(1))
else
    fillfields(handles,handles.defaultjob)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_path_imagedirectory_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_path_imagedirectory_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_imagedirectory_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_imagedirectory_ed as a double
handles = guidata(hObject);


imagedir = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).imagedirectory =  imagedir;


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_path_imagename_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_path_imagename_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_imagename_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_imagename_ed as a double
handles = guidata(hObject);


imagename = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).imagename =  imagename;


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes during object creation, after setting all properties.
function GUI_st_path_firstimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_path_firstimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_firstimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_firstimage_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).firstimage =  str2num(numb);

% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_path_lastimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_path_lastimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_lastimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_lastimage_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).lastimage =  str2num(numb);

% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in GUI_st_iq_set_pb.
function GUI_st_iq_set_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_set_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leveldeterminer(hObject);
handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
fillFields(handles,handles.jobs(projNum))

% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_background_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_iq_fi_background_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_background_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_background_ed as a double
handles = guidata(hObject);



% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_nucleus_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_iq_fi_nucleus_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_nucleus_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_nucleus_ed as a double
handles = guidata(hObject);


% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_background_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_iq_la_background_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_background_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_background_ed as a double
handles = guidata(hObject);


% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_nucleus_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_iq_la_nucleus_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_nucleus_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_nucleus_ed as a double
handles = guidata(hObject);



% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_maxsearch_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsearch_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_bp_maxsearch_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsearch_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_maxsearch_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_maxsearch_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String')

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).maxsearch =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_timestepslide_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_timestepslide_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on
% Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_eo_timestepslide_pm.
function GUI_st_eo_timestepslide_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_timestepslide_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_timestepslide_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_timestepslide_pm
handles = guidata(hObject);



% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_mmpixel_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_mmpixel_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_bp_mmpixel_pm.
function GUI_st_bp_mmpixel_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_mmpixel_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_bp_mmpixel_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_bp_mmpixel_pm
handles = guidata(hObject);



% Update handles structure
guidata(hObject, handles);


%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_bp_minsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_bp_minsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_minsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_minsize_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).minsize =  str2num(numb);



% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_maxsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_bp_maxsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_maxsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_maxsize_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).maxsize =  str2num(numb);



% Update handles structure
guidata(hObject, handles);


%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_minsdist_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_bp_minsdist_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_minsdist_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_minsdist_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).minsdist =  str2num(numb);



% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_setminsize_pb.
function GUI_st_bp_setminsize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setminsize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SetCellValues(hObject,1);
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_setmaxsize_pb.
function GUI_st_bp_setmaxsize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setmaxsize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


SetCellValues(hObject,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_setmindist_pb.
function GUI_st_bp_setmindist_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setmindist_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);


SetCellValues(hObject,3);


% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_minedge_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_minedge_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_eo_minedge_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_minedge_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_minedge_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_minedge_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).minedge =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_eo_noiseparameter_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_noiseparameter_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_eo_noiseparameter_pm.
function GUI_st_eo_noiseparameter_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_noiseparameter_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_noiseparameter_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_noiseparameter_pm
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).noiseparameter =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_mincorrqualtempl_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mincorrqualtempl_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_eo_mincorrqualtempl_pm.
function GUI_st_eo_mincorrqualtempl_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mincorrqualtempl_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_mincorrqualtempl_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_mincorrqualtempl_pm
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).mincorrqualtempl =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_loadsettings_pb.
function GUI_st_bp_loadsettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_loadsettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');


imgDir = handles.jobs(projNum).imagedirectory;
imgNam = handles.jobs(projNum).imagename;
firstImg = handles.jobs(projNum).firstimage;
lastImg = handles.jobs(projNum).lastimage;
increment = handles.jobs(projNum).increment;
savedirectory = handles.jobs(projNum).savedirectory;

cd(handles.jobs(1).savedirectory);

[filename,jobValPath] = uigetfile({'*.mat','mat-files'},'Please select a file named jobvalues.mat');

if ~strcmp(filename,'jobvalues.mat')
   disp('select a file named jobvalues.mat!!!!!!!!!')
   return
end

cd(jobValPath)


load('jobvalues.mat');
handles.jobs(projNum) = jobvalues;

handles.jobs(projNum).imagedirectory = imgDir;
handles.jobs(projNum).imagename = imgNam;
handles.jobs(projNum).firstimage = firstImg;
handles.jobs(projNum).lastimage = lastImg;
handles.jobs(projNum).increment = increment;
handles.jobs(projNum).savedirectory = savedirectory;

fillFields(handles,handles.jobs(projNum))

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_savesettings_pb.
function GUI_st_bp_savesettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_savesettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_bp_defaultsettings_pb.
function GUI_st_bp_defaultsettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_defaultsettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');


imgDir = handles.jobs(projNum).imagedirectory;
imgNam = handles.jobs(projNum).imagename;
firstImg = handles.jobs(projNum).firstimage;
lastImg = handles.jobs(projNum).lastimage;
increment = handles.jobs(projNum).increment;
savedirectory = handles.jobs(projNum).savedirectory;

handles.jobs(projNum) = handles.defaultjob;

handles.jobs(projNum).imagedirectory = imgDir;
handles.jobs(projNum).imagename = imgNam;
handles.jobs(projNum).firstimage = firstImg;
handles.jobs(projNum).lastimage = lastImg;
handles.jobs(projNum).increment = increment;
handles.jobs(projNum).savedirectory = savedirectory;

fillFields(handles,handles.jobs(projNum))

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_leveladjust_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_leveladjust_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_eo_leveladjust_pm.
function GUI_st_eo_leveladjust_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_leveladjust_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_leveladjust_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_leveladjust_pm
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).leveladjust =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in GUI_st_run_pb.
function GUI_st_run_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_run_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

jobList = get(handles.GUI_st_job_lb,'String');
if iscell(jobList)
    howmanyjobs = length(jobList); 
else
    return
end

for projNum = 1:howmanyjobs
        Increment = handles.jobs(projNum).increment;

		possibleImg = handles.jobs(projNum).firstimage;
		while (possibleImg+Increment) <= handles.jobs(projNum).lastimage
              possibleImg = possibleImg+Increment;
		end
		if handles.jobs(projNum).lastimage > possibleImg
		     handles.jobs(projNum).lastimage = possibleImg;
        end

            
		%%%%% save the definite version of jobvalues %%%%
		cd(handles.jobs(projNum).savedirectory)
		jobvalues = handles.jobs(projNum);
		save ('jobvalues','jobvalues')
		clear jobvalues
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 		try
           trackmater(hObject,projNum);
% 		catch    
%            disp(['job number ',num2str(projNum),' had an error and could not be completed'])
%            disp(lasterr)
% 		end
           
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);



% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_eo_sizetemplate_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_sizetemplate_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_eo_sizetemplate_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_sizetemplate_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_sizetemplate_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_sizetemplate_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).sizetemplate =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_mintrackcorrqual_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mintrackcorrqual_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_eo_mintrackcorrqual_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mintrackcorrqual_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_mintrackcorrqual_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_mintrackcorrqual_ed as a double
handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).mintrackcorrqual = str2num(numb);

% Update handles structure
guidata(hObject, handles);


%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function GUI_st_path_savedirectory_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_path_savedirectory_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_savedirectory_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_savedirectory_ed as a double
handles = guidata(hObject);


savedirectory = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).savedirectory =  savedirectory;


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in GUI_st_path_savedirectory_browse_pb.
function GUI_st_path_savedirectory_browse_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_browse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

handles = guidata(hObject);

savedirectory = uigetdir;

if savedirectory== 0
    return
end

set(handles.GUI_st_path_savedirectory_ed,'String',savedirectory);

projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).savedirectory =  savedirectory;



% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_path_increment_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_increment_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function GUI_st_path_increment_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_increment_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_increment_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_increment_ed as a double
handles = guidata(hObject);


numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).increment =  str2num(numb);


% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues = handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton22.
function GUI_st_test_pb_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
testbutton(hObject);

handles = guidata(hObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_halolevel_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_iq_fi_halolevel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_halolevel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_halolevel_ed as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_halolevel_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_iq_la_halolevel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_halolevel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_halolevel_ed as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function GUI_st_bitdepth_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bitdepth_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in GUI_st_bitdepth_pm.
function GUI_st_bitdepth_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bitdepth_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_bitdepth_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_bitdepth_pm
handles = guidata(hObject)


bitDepth = (get(hObject,'Value')*2)+6;

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).intensityMax =  2^bitDepth-1;
handles.jobs(projNum).bitdepth = bitDepth;
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes during object creation, after setting all properties.
function GUI_st_path_timeperframe_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_timeperframe_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_st_path_timeperframe_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_timeperframe_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_timeperframe_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_timeperframe_ed as a double

handles = guidata(hObject);

numb = get(hObject,'String');

%select current project
projNum = get(handles.GUI_st_job_lb,'Value');

handles.jobs(projNum).timeperframe =  str2num(numb);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
