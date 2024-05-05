function varargout = theOtherChannelReadingProcessGUI(varargin)
% THEOTHERCHANNELREADINGPROCESSGUI M-file for theOtherChannelReadingProcessGUI.fig
%      THEOTHERCHANNELREADINGPROCESSGUI, by itself, creates a new THEOTHERCHANNELREADINGPROCESSGUI or raises the existing
%      singleton*.
%
%      H = THEOTHERCHANNELREADINGPROCESSGUI returns the handle to a new THEOTHERCHANNELREADINGPROCESSGUI or the handle to
%      the existing singleton*.
%
%      THEOTHERCHANNELREADINGPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THEOTHERCHANNELREADINGPROCESSGUI.M with the given input arguments.
%
%      THEOTHERCHANNELREADINGPROCESSGUI('Property','Value',...) creates a new THEOTHERCHANNELREADINGPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before theOtherChannelReadingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to theOtherChannelReadingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help theOtherChannelReadingProcessGUI

% Last Modified by GUIDE v2.5 10-Jun-2018 16:44:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @theOtherChannelReadingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @theOtherChannelReadingProcessGUI_OutputFcn, ...
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


% --- Executes just before theOtherChannelReadingProcessGUI is made visible.
function theOtherChannelReadingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',0);

% Parameter setup
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;

% Get list of samplable output
samplableInput = TheOtherChannelReadingProcess.getSamplableInput;
samplableProcessNames = {samplableInput.processName};
samplableOutput = {samplableInput.samplableOutput};

% Read available samplable  output for the given movie
imageProc = userData.MD.processes_;
samplableProcessIndex = cellfun(@(x)userData.MD.getProcessIndex(x,Inf),...
    unique({samplableInput.processName}),'UniformOutput',false);
samplableProcessIndex = sort(vertcat(samplableProcessIndex{:}),'descend');
samplableChannelIndex = arrayfun(@(x)find(imageProc{x}.checkChannelOutput),...
    samplableProcessIndex,'UniformOutput',false);
samplableType = @(pid) cellfun(@(x)isa(imageProc{pid},x),samplableProcessNames);
samplableProcessOutput = arrayfun(@(x) samplableOutput(samplableType(x)),...
    samplableProcessIndex,'UniformOutput',false);

fontSize=6.5;
% Create templates for GUI generation
createText= @(pos,text) uicontrol(handles.uipanel_samplableInput,'Style','text',...
    'Position',[40 pos 350 20],'String',text,'HorizontalAlignment','left','FontSize',fontSize);
createChannelBox= @(i,j,k,pos,varargin) uicontrol(handles.uipanel_samplableInput,'Style','checkbox',...
    'Position',[400+30*k pos 20 20],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j) '_channel' num2str(k)]);

% Create checkboxes for samplable processes
hPosition1 =10;
for i = 1:numel(samplableProcessIndex)
    iProc = samplableProcessIndex(i);
    output=imageProc{iProc}.getDrawableOutput;
    for j=numel(samplableProcessOutput{i}):-1:1
        samplableOutputVar = samplableProcessOutput{i}(j);
        iOutput = find(strcmp(samplableOutputVar,{output.var}));
        createText(hPosition1,[imageProc{iProc}.getName  ' - ' output(iOutput).name]);
        arrayfun(@(x) createChannelBox(iProc,iOutput,x,hPosition1),samplableChannelIndex{i});
        hPosition1=hPosition1+20;
    end
end

% Create checkbox for samplable raw iamges
createText(hPosition1,'Raw images');
arrayfun(@(x) createChannelBox(0,1,x,hPosition1),1:numel(userData.MD.channels_));
hPosition1=hPosition1+20;

% Add cosmetic title for channels
uicontrol(handles.uipanel_samplableInput,'Style','text',...
    'Position',[250 hPosition1 150 20],'String','Channels','FontSize', fontSize);
arrayfun(@(x) uicontrol(handles.uipanel_samplableInput,'Style','text',...
    'Position',[400+30*x hPosition1 20 20],'String',x,'FontSize',fontSize),1:numel(userData.MD.channels_));

% Fix GUI position/size
a=get(get(handles.uipanel_samplableInput,'Children'),'Position');
P=vertcat(a{:});
panelSize = [max(P(:,1)+P(:,3))+10 max(P(:,2)+P(:,4))+20];
pos = get(handles.uipanel_samplableInput,'Position');
dh= panelSize(2) - pos(2);
dL = max(0, panelSize(1) - pos(3));
set(handles.figure1,'Position',get(handles.figure1,'Position')+[0 -dh dL dh]);
set(handles.uipanel_samplableInput,'Position',get(handles.uipanel_samplableInput,'Position')+[0 0 dL dh]);

% Move various graphic elements
set(handles.axes_help,'Position',get(handles.axes_help,'Position')+[dL dh 0 0]);
set(handles.text_processName,'Position',get(handles.text_processName,'Position')+[0 dh 0 0]);
set(handles.text_copyright,'Position',get(handles.text_copyright,'Position')+[0 dh 0 0]);

% Update handles structure and attach it to the main figure
handles = guihandles(handles.figure1);
guidata(handles.figure1, handles); 

% Read the default sampled processes
initProcessIndex =funParams.ProcessIndex;
initChannelIndex =funParams.ChannelIndex;
initOutputName =funParams.OutputName;

if ~iscell(initProcessIndex),initProcessIndex={initProcessIndex}; end
if ~iscell(initChannelIndex),initChannelIndex={initChannelIndex};end
if ~iscell(initOutputName),initOutputName={initOutputName};end

for i=1:numel(initProcessIndex)
    procId=initProcessIndex{i};
    % Retrieve the id, process nr and channel nr of the selected imageProc
    if ~isempty(procId)
        output=imageProc{procId}.getDrawableOutput;
        iOutput=find(strcmp(initOutputName{i},{output.var}));
        procTag=['checkbox_process' num2str(procId) '_output' num2str(iOutput)];
    else
        procTag='checkbox_process0_output1';
    end
        
    % Set checkbox values to 1
    for j=initChannelIndex{i}(:)'      
        set(handles.([procTag '_channel' num2str(j)]),'Value',1);
    end
end

% % Set up dark current paths
% if ~all(cellfun(@isempty,userData.crtProc.inFilePaths_(:)))
%     set(handles.listbox_darkCurrentChannels, 'String', ...
%         userData.crtProc.inFilePaths_(1,funParams.iChanSlave ));
% end
% set(handles.listbox_darkCurrentChannels, 'Value',1);
% set(handles.checkbox_FAregistration,'Value',funParams.doFAregistration);

userData.userDir=userData.MD.outputDirectory_;

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% uicontrol(handles.pushbutton_done)


% --- Outputs from this function are returned to the command line.
function varargout = theOtherChannelReadingProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Check user input
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

% Retrieve list of checked boxes
h=findobj(handles.figure1,'Style','Checkbox','-and','Value',1,'-not',...
    'Tag','checkbox_applytoall');
tags= get(h,'Tag');
if ischar(tags), tags={tags}; end
tokens=cellfun(@(x)regexp(x,'^checkbox_process(\d+)_output(\d+)_channel(\d+)','tokens'),...
    tags);

% Read process index, channel index and output index
procID=cellfun(@(x)  str2double(x{1}),tokens);
outputID=cellfun(@(x)  str2double(x{2}),tokens);
chanID=cellfun(@(x)  str2double(x{3}),tokens);
uProcID=unique(procID);

% Initialize funParams
funParams.ProcessIndex={};
funParams.ChannelIndex={};
funParams.OutputName={};

% If raw images are selected
if ismember(0,uProcID)
    funParams.ProcessIndex{1}=[];
    funParams.ChannelIndex{1}=sort(chanID(procID==0));
    funParams.OutputName{1}='';
end
% List output to sample
for pid=uProcID(uProcID>0)'
    output = userData.MD.processes_{pid}.getDrawableOutput;
    uOutputID=unique(outputID(procID==pid));
    for j=uOutputID'
        funParams.ProcessIndex{end+1}=pid;
        funParams.ChannelIndex{end+1}=sort(chanID(procID==pid & outputID==j));
        funParams.OutputName{end+1}=output(j).var;
    end
end

% Apply new channel index and new shade images for sanity check
% funParams.iChanMaster = channelIndex; 
% funParams.iChanSlave = iSlaveChannels; 
funParams.doFAregistration = get(handles.checkbox_FAregistration,'Value');

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% % Call back function of 'Apply' button
% userData = get(handles.figure1, 'UserData');
% if isempty(userData), userData = struct(); end
% channelProps = get (handles.listbox_selectedChannels, {'Userdata','String'});
% channelIndex=channelProps{1};
% channelString=channelProps{2};
% slaveChannels = get(handles.listbox_darkCurrentChannels, 'String');
% % Check user input
% if isempty(slaveChannels)
%     errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
%     return;
% elseif numel(slaveChannels)>1
%     errordlg('Please select only one input channel from ''Available Channels''.','Setting Error','modal')
%     return;
% end
% iSlaveChannels=find(strcmp(userData.crtProc.inFilePaths_,slaveChannels));
% 
% 
% % Save old process settings (for recovery if sanity check fails)
% funParams = userData.crtProc.funParams_;
% 
% % Apply new channel index and new shade images for sanity check
% funParams.iChanMaster = channelIndex; 
% funParams.iChanSlave = iSlaveChannels; 
% funParams.doFAregistration = get(handles.checkbox_FAregistration,'Value');
% 
% parseProcessParams(userData.crtProc,funParams);
% 
% % Set parameters
% processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
% 
% % --- Executes on button press in pushbutton_deleteDarkCurrentChannel.
% function pushbutton_deleteDarkCurrentChannel_Callback(hObject, eventdata, handles)
% % Call back function of 'delete' button
% contents = get(handles.listbox_darkCurrentChannels,'String');
% % Return if list is empty
% if isempty(contents)
%     return;
% end
% num = get(handles.listbox_darkCurrentChannels,'Value');
% 
% % Delete selected item
% contents(num) = [ ];
% 
% % Refresh listbox
% set(handles.listbox_darkCurrentChannels,'String',contents);
% % Point 'Value' to the second last item in the list once the 
% % last item has been deleted
% if (num>length(contents) && num>1)
%     set(handles.listbox_darkCurrentChannels,'Value',length(contents));
% end

guidata(hObject, handles);


% --- Executes on button press in checkbox_gaussianfilter.
function checkbox_gaussianfilter_Callback(hObject, eventdata, handles)

if get(hObject, 'Value'), state='on'; else state='off'; end
set([handles.text_filters1 handles.edit_sigma], 'Enable',state);

% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)

% call back of 'Up' button
id = get(handles.listbox_darkCurrentChannels,'Value');
contents = get(handles.listbox_darkCurrentChannels,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

set(handles.listbox_darkCurrentChannels, 'string', contents);
set(handles.listbox_darkCurrentChannels, 'value', id-1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in checkbox_FAregistration.
function checkbox_FAregistration_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_FAregistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_FAregistration
