function varargout = fsmGuiMain(varargin)
% fsmGuiMain Application M-file for fsmGuiMain.fig
%    FIG = fsmGuiMain launch fsmGuiMain GUI.
%    fsmGuiMain('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 11-Mar-2005 16:45:29

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  INITIALIZING WINDOW
%  LOADING DEFAULT VALUES (fsmParams), SETTING GUI FIELDS WITH DEFAULT VALUES
%  STORING DEFAULT VALUES INTO defaultFsmParams
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0  % LAUNCH GUI

    % Check whether fsmGuiMain is already running
    if ~isempty(findall(0,'Tag','fsmGuiMain','Name','SpeckTackle'))
        alreadyOpen=1;
    else
        alreadyOpen=0;
    end
    
	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % INITIALIZE USER INTERFACE IF IT WAS NOT OPEN BEFORE
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Fix the weird background color in Linux
    if isunix
        set(handles.fsmGuiMain,'Color',[0.701961 0.701961 0.701961]);
        set(handles.numberEdit,'BackgroundColor',[0.701961 0.701961 0.701961]);
        set(handles.bitDepthEdit,'BackgroundColor',[0.701961 0.701961 0.701961]);
        %set(handles.textPsfSigma,'BackgroundColor',[0.701961 0.701961 0.701961]);
        set(handles.editThreshold,'BackgroundColor',[0.701961 0.701961 0.701961]);
        set(handles.expPopup,'BackgroundColor',[0.701961 0.701961 0.701961]);        
    end
    
    if alreadyOpen==0
        
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
        if exist(fsmParamPath, 'file') == 2
            load(fsmParamPath); % Load defaut values stored in the fsm directory
        else
            % Initialize fsmParam with default values
            fsmParam=fsmGetParamDflts;
        end
        
        %Since we are going to distinguish the various sigma in fsm from now 
        %(Mar. 11, 2005) on, for backward compatiblity, I copy the old 'sigma'
        %to 'filterSigma'.
        if ~isfield(fsmParam.prep,'filterSigma')
           fsmParam.prep.filterSigma = fsmParam.prep.sigma;
        end
        if ~isfield(fsmParam.prep,'psfSigma')
           fsmParam.prep.psfSigma = fsmParam.prep.sigma;
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
            if ~exist(fsmExpParamPath, 'file')
                % No database found in user-defined directory
                % Reverting to default database
                fsmExpParamPath=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
            end
        end
        if ~exist(fsmExpParamPath, 'file')
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
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CHECK WHETHER A PROJECT IS OPEN IN FSMCENTER
    %
    %    Anyway, open fsmCenter if it not yet open
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update if needed
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
    if ~isempty(projDir)
        % Set up the project path
        if strcmp(projDir(end),filesep)==1
            workPath=[projDir,char(settings.subProjects(1)),filesep];
        else
            workPath=[projDir,filesep,char(settings.subProjects(1)),filesep];
        end
        if iscell(settings.imageDir)
            imagePath=settings.imageDir{1};
        else
            imagePath=settings.imageDir;
        end
        
        % Store project information if fsmGuiMain's UserData
        set(handles.fsmGuiMain,'UserData',settings);
        
    else
        workPath=[];
        imagePath=[];
    end
    set(handles.pathEdit,'String',workPath);
    set(handles.textImage,'String',imagePath);

    fsmParam = get(handles.start,'UserData');
    %Get physical parameters from fsmCenter and broadcast them to 'fsmParam'.
    physiParam = handlesFsmCenter.physiParam;
    fsmParam.prep.psfSigma    = physiParam.psfSigma;
    fsmParam.prep.filterSigma = physiParam.psfSigma;

    set(handles.start,'UserData',fsmParam);
    set(handles.defaultButton,'UserData',fsmParam);
    fsmGuiWriteParameters(fsmParam,handles);

    % Update the GUI if a fsmParam.mat file exists in the workPath
    catchPathChange(workPath,handles,settings,fsmParam);
    
    %Check whether the image path stored in fsmParam is compatible with the
    %one in 'settings' which is setup in fsmCenter. The incompatibility can
    %be caused by platform dependent directory format.
    fsmParam = get(handles.start,'UserData');
    if iscell(settings.imageDir)
        imageDir = settings.imageDir{1};
    else
        imageDir = settings.imageDir;
    end
    noProblem = 1;
    if isempty(fsmParam.main.imagePath)
        fsmParam.main.imagePath = imageDir;
    end
    if ~samdir(fsmParam.main.imagePath,imageDir)
        if ~isdir(fsmParam.main.imagePath)
            %This could be a platform problem.
            if ispc == 1
                if ~isempty(fsmParam.main.imagePath) && strcmp(fsmParam.main.imagePath(1),'/')
                    %This is a Unix directory. Try to convert it to PC
                    %format.
                    imgDrive = getDriveName(imageDir);
                    imagePath = dirUnix2PC(fsmParam.main.imagePath,imgDrive);
                    if ~samdir(imagePath,imageDir)
                        %Still not same dir. Image path has been changed
                        %since last tracking. Give warning.
                        noProblem = 0;
                    else
                        fsmParam.main.imagePath = imagePath;
                    end
                else
                    noProblem = 0;
                end
            elseif isunix == 1
                if ~isempty(fsmParam.main.imagePath) && ~strcmp(fsmParam.main.imagePath(1),'/')
                    imgDrive = getDriveName(imageDir);
                    imagePath = dirPC2Unix(fsmParam.main.imagePath,imgDrive);
                    if ~samdir(imagePath,imageDir)
                        %Still not same dir. Image path has been changed
                        %since last tracking. Give warning.
                        noProblem = 0;
                    else
                        fsmParam.main.imagePath = imagePath;
                    end
                else
                    noProblem = 0;                    
                end
            end
        else
            %It is not platform problem. Image path has been changed since
            %last tracking.
            noProblem = 0;
        end
    end
    
    if ~noProblem
        %Ask the user to make a choice.
        question = sprintf('%s\n%s\n%s\n%s\n%s\n%s','The image path from project setup: ', ...
            imageDir, ...
            'is different from the stored image path of last specktackle tracking:', ...
            [fsmParam.main.imagePath '.'], 'Do you want to continue with the project image path:', ...
            [imageDir '?']);
        answer = questdlg(question,'Alert','Yes','No','Yes');
        
        if strcmp(answer,'Yes') == 1
            fsmParam.main.imagePath = imageDir;
        else
            fsmH=findall(0,'Tag','fsmGuiMain','Name','SpeckTackle'); % Get the handle of specktackle.
            if ishandle(fsmH)
                delete(fsmH);
                return;
            end
        end
    end
    set(handles.start,'UserData',fsmParam);
    fsmGuiWriteParameters(fsmParam,handles);

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

% Check the selection of modules
sel=[get(handles.checkPrepModule,'Value'), ...
        get(handles.checkTrackModule,'Value'), ...
        get(handles.checkBuildModule,'Value'), ...
        get(handles.checkKinModule,'Value'), ...
        get(handles.checkDispModule,'Value')];
indx0=find(sel==0);
indx1=find(sel==1);
if any([indx0>min(indx1) & indx0<max(indx1)])
    uiwait(msgbox('Please check the selection of your modules.','Error','Error','modal'));
    return
end    

% Read fsmParam from start button
fsmParam=get(handles.start,'UserData');

% Read project structure from fsmGuiMain's UserData and sotre it in fsmParam.project
settings=get(handles.fsmGuiMain,'UserData');
if isempty(settings)
    error('Could not retrieve project settings from fsmGuiMain''s ''UserData'' field.');
end
projectDir=settings.projDir;
subProjects=settings.subProjects;
fsmParam.project.path=projectDir;
fsmParam.project.corr=subProjects{strmatch('corr',subProjects)};
fsmParam.project.edge=subProjects{strmatch('edge',subProjects)};
fsmParam.project.lpla=subProjects{strmatch('lpla',subProjects)};
fsmParam.project.merg=subProjects{strmatch('merg',subProjects)};
fsmParam.project.post=subProjects{strmatch('post',subProjects)};
fsmParam.project.tack=subProjects{strmatch('tack',subProjects)};

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
% AUTO POLYGON CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = autoPolCheck_Callback(h, eventdata, handles, varargin)
% Nothing has to be done

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PATH SELECTION CALLBACKS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function varargout = ButtonBrowse_Callback(h, eventdata, handles, varargin)
% 
% % Select path
% newpath=uigetdir('','Select work directory');
%
% % Check whether an fsmParam.mat file already exists
% catchPathChange(newpath,handles);
% 
% %%%%%
% 
% function pathEdit_Callback(hObject, eventdata, handles)
% 
% newpath=get(handles.pathEdit,'String');
% % Check whether an fsmParam.mat file already exists
% catchPathChange(newpath,handles);


%%%%%

function catchPathChange(newpath,handles,settings,inFsmParam)

if newpath~=0

    % If an fsmParam.mat file is found in this path, use this as fsmParam
    if exist([newpath,filesep,'fsmParam.mat'], 'file') == 2
        
        % Inform user
        uiwait(msgbox('A file ''fsmParam.mat'' has been found in the selected path. Loading and setting stored values.','Info','help','modal'));
        
        % Load the file found
        load([newpath,filesep,'fsmParam.mat']);  

        %For backward compatibility.
        if ~isfield(fsmParam.prep,'filterSigma')
           fsmParam.prep.filterSigma = fsmParam.prep.sigma;
        end
        if ~isfield(fsmParam.prep,'psfSigma') && nargin == 4 && ...
           isfield(inFsmParam.prep,'psfSigma')
           fsmParam.prep.psfSigma = inFsmParam.prep.psfSigma;
        end
            
    else
        
        % Get default values
        fsmParam=get(handles.defaultButton,'UserData');
        
    end

    % Fill all fields with the values from fsmParam
    fsmGuiWriteParameters(fsmParam,handles);
    
    % Link fsmParam to the start button (the default values are not modified)
    set(handles.start,'UserData',fsmParam);
    
    if isempty(fsmParam.main.path)
        
        % Update user interface with project paths
        set(handles.pathEdit,'String',newpath);
        set(handles.textImage,'String',char(settings.imageDir));

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREPROCESSING MODULE ON/OFF CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = checkPrepModule_Callback(h, eventdata, handles, varargin)

if get(handles.checkPrepModule,'Value')==0
    togglePrepModule(handles,0);
else
    togglePrepModule(handles,1);
end
    
function togglePrepModule(handles,value)

if value==0
    set(handles.autoPolCheck,'Enable','off');
    set(handles.textCameraCalPar,'Enable','off');
    set(handles.expPopup,'Enable','off');
    set(handles.textAdvancedPrep,'Enable','off');
    set(handles.textGauss,'Enable','off');
    set(handles.editGauss,'Enable','off');
    set(handles.subpixel,'Enable','off');
    set(handles.primaryRadio,'Enable','off');
    set(handles.tertiaryRadio,'Enable','off');
    set(handles.drawROICheck,'Enable','off');
    set(handles.loadROICheck,'Enable','off');   
    set(handles.textDescr,'Enable','off');
    set(handles.textSigma,'Enable','off');
    set(handles.editSigma,'Enable','off');
    set(handles.percText,'Enable','off');
    set(handles.percEdit,'Enable','off');
    set(handles.orderText,'Enable','off');
    set(handles.orderEdit,'Enable','off');
else
    set(handles.autoPolCheck,'Enable','on');
    set(handles.textCameraCalPar,'Enable','on');
    set(handles.expPopup,'Enable','on');
    set(handles.textAdvancedPrep,'Enable','on');
    set(handles.textGauss,'Enable','on');
    set(handles.subpixel,'Enable','on');
    set(handles.editGauss,'Enable','on');
    set(handles.primaryRadio,'Enable','on');
    set(handles.tertiaryRadio,'Enable','on');
    set(handles.drawROICheck,'Enable','on');
    set(handles.loadROICheck,'Enable','on');
    set(handles.textDescr,'Enable','on');
    set(handles.textSigma,'Enable','on');
    set(handles.editSigma,'Enable','on');    
    set(handles.percText,'Enable','on');
    set(handles.percEdit,'Enable','on');
    set(handles.orderText,'Enable','on');
    set(handles.orderEdit,'Enable','on');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% KINETIC ANALYSIS MODULE ON/OFF CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = checkKinModule_Callback(h, eventdata, handles, varargin)

if get(handles.checkKinModule,'Value')==0
    toggleKinModule(handles,0);
else
    toggleKinModule(handles,1);
end

function toggleKinModule(handles,value)

if value==0
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

if get(handles.checkTrackModule,'Value')==0
    fsmGuiMain('toggleTrackModule',handles,0);
else
    fsmGuiMain('toggleTrackModule',handles,1);
    
    % Update fields depending on which tracker is chosen
    fsmGuiMain('popupTracker_Callback',hObject, eventdata, handles);
    
end    

% This function turns everything on|off in the tracking module
function toggleTrackModule(handles,value)

if value==0
    set(handles.textThreshold,'Enable','off');
    set(handles.editThreshold,'Enable','off');
    set(handles.checkEnhTrack,'Enable','off');
    set(handles.checkGrid,'Enable','off');
    set(handles.textInfluence,'Enable','off');
    set(handles.editInfluence,'Enable','off');
    set(handles.textCorrLength,'Enable','off');
    set(handles.editCorrLength,'Enable','off');
    set(handles.popupTracker,'Enable','off');
    set(handles.textTracker,'Enable','off');
    set(handles.popupTrackInit,'Enable','off');
    set(handles.checkTrackInit,'Enable','off');
    set(handles.textInitRadius,'Enable','off');
    set(handles.editInitRadius,'Enable','off');    
else
    set(handles.textThreshold,'Enable','on');
    set(handles.editThreshold,'Enable','on');
    set(handles.checkEnhTrack,'Enable','on');
    set(handles.checkGrid,'Enable','on');
    set(handles.textInfluence,'Enable','on');
    set(handles.editInfluence,'Enable','on');
    set(handles.textCorrLength,'Enable','on');
    set(handles.editCorrLength,'Enable','on');
    set(handles.popupTracker,'Enable','on');
    set(handles.textTracker,'Enable','on');
    set(handles.popupTrackInit,'Enable','on');
    set(handles.checkTrackInit,'Enable','on');
    set(handles.textInitRadius,'Enable','on');
    set(handles.editInitRadius,'Enable','on');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SAVE PARAMETERS (FOR BATCH JOBS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveParamsButton_Callback(hObject, eventdata, handles)
userPath=get(handles.pathEdit,'String');
if isempty(userPath)
    uiwait(msgbox('Please setup a project in fsmCenter.','Warning','warn'));
    return;
else
    % Create ROOT directory
    if ~exist(userPath, 'dir')
        % Directory does not exist - create it
        if userPath(2)==':';   % Windows
            % Drive letter specified
            mkdir(userPath(1:3),userPath(4:end));
        else
            mkdir(userPath);
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
fsmGuiUpdateConfidences(1);


% --------------------------------------------------------------------
function varargout = primaryRadio_Callback(h, eventdata, handles, varargin)

set(handles.primaryRadio,'Value',1);
set(handles.tertiaryRadio,'Value',0);
set(handles.percText,'Enable','off');
set(handles.percEdit,'Enable','off');
set(handles.orderText,'Enable','off');
set(handles.orderEdit,'Enable','off');

% --------------------------------------------------------------------
function varargout = tertiaryRadio_Callback(h, eventdata, handles, varargin)

set(handles.primaryRadio,'Value',0);
set(handles.tertiaryRadio,'Value',1);
set(handles.percText,'Enable','on');
set(handles.percEdit,'Enable','on');
set(handles.orderText,'Enable','on');
set(handles.orderEdit,'Enable','on');

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

% Since we need to change 'filterSigma', we also get 'fsmParam'.
fsmParam = get(handles.start,'UserData');

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

% Write gauss kernel
if exp==1
    set(handles.editSigma,'String',num2str(defaultFsmParam.prep.sigma));
else
    if ~isfield(fsmExpParam(exp-1),'caliSigma') | isempty(fsmExpParam(exp-1).caliSigma)
        set(handles.editSigma,'String',num2str(defaultFsmParam.prep.sigma));
    else
        set(handles.editSigma,'String',num2str(fsmExpParam(exp-1).caliSigma));
    end
end

% Write description
if exp==1
    set(handles.textDescr,'String','Experiment description');
else
    set(handles.textDescr,'String',fsmExpParam(exp-1).description);
end

%Update 'filterSigma' in 'fsmParam'.
filterSigma = str2double(get(handles.editSigma,'String'));
fsmParam.prep.filterSigma = filterSigma;
set(handles.start,'UserData',fsmParam);

% Check whether it is an optimized experiment or not
if exp~=1
    if fsmExpParam(exp-1).quantile==0
        
        % Not optimized
        if strcmp(get(handles.editZValue,'Enable'),'off')
            set(handles.editZValue,'String','1.96');
            fsmGuiUpdateConfidences(1);
        end
    else
        % Optimized
        set(handles.editZValue,'String',num2str(fsmExpParam(exp-1).quantile));
        fsmGuiUpdateConfidences(0);
    end
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
if get(handles.drawROICheck,'Value')==1
    set(handles.loadROICheck,'Value',0);
end

% --------------------------------------------------------------------
function fsmGuiMain_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fsmH=findall(0,'Tag','fsmGuiMain','Name','SpeckTackle'); % Get the handle of fsmPostProc

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
if fsmH~=hObject && force==1
    delete(fsmH); % Force closing
else
    choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
    switch choice,
        case 'Yes', delete(fsmH);
        case 'No', return;
    end
end


% --- Executes during object deletion, before destroying properties.
function fsmGuiMain_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fsmGuiMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushCalReload.
function pushCalReload_Callback(hObject, eventdata, handles)
% hObject    handle to pushCalReload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function editInfluence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editInfluence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function editInfluence_Callback(hObject, eventdata, handles)
% hObject    handle to editInfluence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editInfluence as text
%        str2double(get(hObject,'String')) returns contents of editInfluence as a double


% --- Executes during object creation, after setting all properties.
function editSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function editSigma_Callback(hObject, eventdata, handles)
% hObject    handle to editSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSigma as text
%        str2double(get(hObject,'String')) returns contents of editSigma as a double

fsmParam = get(handles.start,'UserData');
filterSigma = str2double(get(hObject,'String'));
if isempty(filterSigma)
   set(handles.editSigma,'String',num2str(fsmParam.prep.filterSigma));
   errordlg('Not valid numerical value','Error','modal');
   return;
end

fsmParam.prep.filterSigma = filterSigma;
set(handles.start,'UserData',fsmParam);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editCorrLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCorrLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function editCorrLength_Callback(hObject, eventdata, handles)
% hObject    handle to editCorrLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCorrLength as text
%        str2double(get(hObject,'String')) returns contents of editCorrLength as a double


% --- Executes on button press in checkTrackInit.
function checkTrackInit_Callback(hObject, eventdata, handles)
% hObject    handle to checkTrackInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkTrackInit
if get(handles.checkTrackInit,'Value')==1
    set(handles.popupTrackInit,'Enable','on');
    fsmGuiMain('popupTrackInit_Callback',handles.popupTrackInit,[],handles);
else
    set(handles.popupTrackInit,'Enable','off');
    fsmGuiMain('popupTrackInit_Callback',handles.popupTrackInit,[],handles);
end

% --- Executes during object creation, after setting all properties.
function popupTracker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupTracker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupTracker.
function popupTracker_Callback(hObject, eventdata, handles)
% hObject    handle to popupTracker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupTracker contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupTracker
value=get(handles.popupTracker,'Value');
if value==1
    % that is the Neural Network Tracker
    fsmGuiMain('toggleTrackModule',handles,1);
    fsmGuiMain('checkEnhTrack_Callback',handles.checkEnhTrack,[],handles); % Turns on|off gird depending on 'iterative'
    fsmGuiMain('checkTrackInit_Callback',handles.checkTrackInit,[],handles); % Turns on|off Initializer
end
if value==2
    % that is Ge's Linear Assignment Tracker Lap
    % matthias changes October
    fsmGuiMain('toggleTrackModule',handles,0);
    set(handles.checkEnhTrack,'Enable','on');
    set(handles.textCorrLength,'Enable','on');
    set(handles.editCorrLength,'Enable','on');
    set(handles.textThreshold,'Enable','on');
    set(handles.editThreshold,'Enable','on');
    set(handles.popupTracker,'Enable','on');
    set(handles.checkTrackInit,'Enable','on');
    set(handles.textTracker,'Enable','on');
    set(handles.textThreshold,'Enable','on');
    set(handles.editThreshold,'Enable','on');   
    set(handles.checkGrid,'Enable','on');   
    fsmGuiMain('checkTrackInit_Callback',handles.popupTrackInit,[],handles);
    fsmGuiMain('popupTrackInit_Callback',handles.popupTrackInit,[],handles);    
end
if value==3
    fsmGuiMain('toggleTrackModule',handles,0);
    set(handles.textTracker,'Enable','on');
    set(handles.popupTracker,'Enable','on');
    set(handles.popupTrackInit,'Value',2);
    set(handles.checkTrackInit,'Enable','on');
    set(handles.checkTrackInit,'Value',1);
    fsmGuiMain('checkEnhTrack_Callback',handles.checkEnhTrack,[],handles); % Turns on|off gird depending on 'iterative'
    fsmGuiMain('checkTrackInit_Callback',handles.checkTrackInit,[],handles); % Turns on|off Initializer
    fsmGuiMain('popupTrackInit_Callback',handles.popupTrackInit,[],handles);
end   

% --- Executes on button press in checkInitKymo.
function checkInitKymo_Callback(hObject, eventdata, handles)
% hObject    handle to checkInitKymo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkInitKymo

% --- Executes during object creation, after setting all properties.
function popupTrackInit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupTrackInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupTrackInit.
function popupTrackInit_Callback(hObject, eventdata, handles)
% hObject    handle to popupTrackInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupTrackInit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupTrackInit
value=get(handles.popupTrackInit,'Value');
if value==1
    % correlation
%     if get(handles.popupTracker,'Value')==2
%         uiwait(warndlg('The Linear Programming tracker does not support initialization by ''Correlation''.','Warning','modal'));
%         set(handles.popupTrackInit,'Value',2)
%         return;
%     end
    if get(handles.popupTracker,'Value')==3
        set(handles.popupTracker,'Value',1)
        fsmGuiMain('popupTracker_Callback',handles.checkTrackInit,[],handles)
    end
    set(handles.textInitRadius,'Enable','off');
    set(handles.editInitRadius,'Enable','off');    
elseif value==2
    % TFT 
    if get(handles.checkTrackInit,'Value')==1
        set(handles.textInitRadius,'Enable','on');
        set(handles.editInitRadius,'Enable','on');        
    else
        set(handles.textInitRadius,'Enable','off');
        set(handles.editInitRadius,'Enable','off');                
    end
else
    error('Impossible value');
end


function loadROICheck_Callback(hObject, eventdata, handles)
if get(handles.loadROICheck,'Value')==1
    set(handles.drawROICheck,'Value',0);
end

function editInitRadius_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editInitRadius_Callback(hObject, eventdata, handles)

function checkTest_Callback(hObject, eventdata, handles)

function subpixel_Callback(hObject, eventdata, handles)

function textPsfSigma_CreateFcn(hObject, eventdata, handles)



