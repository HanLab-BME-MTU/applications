function varargout = spotfinder(varargin)
% SPOTFINDER Application M-file for spotfinder.fig
%    FIG = SPOTFINDER launch spotfinder GUI.
%    SPOTFINDER('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 19-Jun-2001 16:34:50

if nargin == 0  % LAUNCH GUI

    load_sp_constants;
	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    %create status struct
    status.fname='';               % movie filename
    status.datapath=tempdir;       % result path
    status.filter=0;               % filter status (0=not filtered)
    status.spots=0;                % spots status  (0=spots not detected)
    status.multispots=0;           % mulitplespots status (0=spots not detected )
    status.classification=0;       % classifaction status (0= not classified) 
    setappdata(fig,'stat',status);
  
	if nargout > 0
		varargout{1} = fig;
	end

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

% --------------------------------------------------------------------
function varargout = Loadr3d_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the load menu.
status=getappdata(handles.spotfindwin,'stat');
logs=get(handles.logwin,'String');
[movie, fname]= r3dread;
if ~isempty(fname)
    setappdata(handles.spotfindwin,'movie',movie);
    status.fname=fname;
    
    %select save folder
    [filename, dpath] = uiputfile([fname '_'], 'Select folder to save data');
    if(dpath(1)~=0)
        st=mkdir(dpath,[filename '_data']);
        status.datapath=[dpath filename '_data' filesep];
    end;
    
    % reset status
    status.filter=0;
    status.spots=0;
    status.classification=0;
    status.multispots=0;
    logs=['File: ' fname ' loaded'];
else
    logs=strvcat(logs,['No file loaded']);
end;
set(handles.logwin,'String',logs);

%save log
fid = fopen([status.datapath 'status.log'],'w');
clog=cellstr(logs);
fprintf(fid,'%s\n',clog{:});
fclose(fid);

setappdata(handles.spotfindwin,'stat',status);

% --------------------------------------------------------------------
function varargout = LoadData_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the Load data menu.

status=getappdata(handles.spotfindwin,'stat');
[filename, dpath] = uigetfile('', 'Select data file');
if(dpath(1)~=0)
    load([dpath filename]);
    logs=textread([dpath 'status.log'],'%s','whitespace','\n');
    logs=char(logs);
    status.datapath=dpath;
    if exist('fmov')
        %store filtered movie
        setappdata(handles.spotfindwin,'movie',fmov);
        status.filter=1;
    end;
    if exist('cord')
        setappdata(handles.spotfindwin,'clist',cord);
        status.multispots=1; 
        status.spots=1;
    end;
    if exist('slist')
        setappdata(handles.spotfindwin,'slist',slist);
        status.classification=1;
    end;
    set(handles.logwin,'String',logs);
    setappdata(handles.spotfindwin,'stat',status);
end;

% --------------------------------------------------------------------
function varargout = filtmovie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.

%CONST DEFINITIONS
global  FILTERPRM;

status=getappdata(handles.spotfindwin,'stat');
movie=getappdata(handles.spotfindwin,'movie');
logs=get(handles.logwin,'String');
if ~isempty(movie)
    if ~status.filter
        % movie filtered
        status.filter=1;
        fmov=filtermovie(movie,FILTERPRM);
        %remove movie to save memory
        rmappdata(handles.spotfindwin,'movie');
        mf=mean(fmov(:));
        % take bg away
        fmov = fmov-mf;
        %store filtered movie
        setappdata(handles.spotfindwin,'movie',fmov);
        logs=strvcat(logs,['Movie filtered']);
        % save filtered movie
        save([status.datapath 'moviedat'],'fmov');
    else
        logs=strvcat(logs,['Movie is already filtered']);
    end;
else
    logs=strvcat(logs,['No movie loaded']);
    errordlg('No movie loaded');
end;
set(handles.logwin,'String',logs);
setappdata(handles.spotfindwin,'stat',status);

%save log
fid = fopen([status.datapath 'status.log'],'w');
clog=cellstr(logs);
fprintf(fid,'%s\n',clog{:});
fclose(fid);


% --------------------------------------------------------------------
function varargout = FindSpots_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the find spots menu.
status=getappdata(handles.spotfindwin,'stat');
logs=get(handles.logwin,'String');
movie=getappdata(handles.spotfindwin,'movie');
if~isempty(movie)
    if(~status.filter)
        ButtonName=questdlg('Movie not filtered, try to find spots anyway?','Warning','Yes','No','No');
        if strcmp(ButtonName,'No')
            return;
        end;
    end;
    if ~status.spots
        [cord,mnp] = spotfind(movie);
        setappdata(handles.spotfindwin,'clist',cord);
        %set log
        logs=strvcat(logs,['Spots detected']);
        status.spots=1;
    else
        logs=strvcat(logs,['Spots already detected']);
    end;
else
    logs=strvcat(logs,['No movie loaded']);
    errordlg('No movie loaded');
end;
set(handles.logwin,'String',logs);    
setappdata(handles.spotfindwin,'stat',status);

%save log
fid = fopen([status.datapath 'status.log'],'w');
clog=cellstr(logs);
fprintf(fid,'%s\n',clog{:});
fclose(fid);

% --------------------------------------------------------------------
function varargout = FindOverlap_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
status=getappdata(handles.spotfindwin,'stat');
logs=get(handles.logwin,'String');
cord=getappdata(handles.spotfindwin,'clist');
movie=getappdata(handles.spotfindwin,'movie');

if~isempty(cord)
    if ~status.multispots
        cord=findoverlap(movie,cord);
        setappdata(handles.spotfindwin,'clist',cord);
        %set log
        logs=strvcat(logs,['Multiple spots detected']);
        status.multispots=1;
        % save spot list
        save([status.datapath 'moviedat'],'cord','-append');
    else
        logs=strvcat(logs,['multiple spots already detected']);
    end;
else
    logs=strvcat(logs,['No spots list']);
    errordlg('No spots list');
end;
set(handles.logwin,'String',logs);
setappdata(handles.spotfindwin,'stat',status);

%save log
fid = fopen([status.datapath 'status.log'],'w');
clog=cellstr(logs);
fprintf(fid,'%s\n',clog{:});
fclose(fid);

% --------------------------------------------------------------------
function varargout = ClassifySpots_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
status=getappdata(handles.spotfindwin,'stat');
logs=get(handles.logwin,'String');
cord=getappdata(handles.spotfindwin,'clist');
if~isempty(cord)
    if ~status.classification
        slist =classifyspot(cord);
        setappdata(handles.spotfindwin,'slist',slist);
        %set log
        logs=strvcat(logs,['Spots classified']);
        status.classification=1;
        % save spot list
        save([status.datapath 'moviedat'],'slist','-append');
        %export spots
        fstat= savespots([status.datapath 'spots'],slist);
    else
        logs=strvcat(logs,['Spots already classified']);
    end;
else
    logs=strvcat(logs,['No spots list']);
    errordlg('No spots list');
end;
set(handles.logwin,'String',logs);
setappdata(handles.spotfindwin,'stat',status);

%save log
fid = fopen([status.datapath 'status.log'],'w');
clog=cellstr(logs);
fprintf(fid,'%s\n',clog{:});
fclose(fid);


% --------------------------------------------------------------------
function varargout = ShowMovie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
movie=getappdata(handles.spotfindwin,'movie');
slist=getappdata(handles.spotfindwin,'slist');
if~isempty(movie)
    if ~isempty(slist)
        r3dview(movie,slist);
    else
        r3dview(movie);
    end;
else
    errordlg('No movie loaded');
end;


% --------------------------------------------------------------------
function varargout = ShowSpots3D_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
movie=getappdata(handles.spotfindwin,'movie');
slist=getappdata(handles.spotfindwin,'slist');
if ~isempty(slist)
    spotview(slist,movie);
else
    errordlg('Spots not classified');
end;


% --------------------------------------------------------------------
function varargout = autoproc_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.listbox1.

%check if loaded

filtmovie_Callback(h, eventdata, handles, varargin);
FindSpots_Callback(h, eventdata, handles, varargin);
FindOverlap_Callback(h, eventdata, handles, varargin);
ClassifySpots_Callback(h, eventdata, handles, varargin);

% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.listbox1.
