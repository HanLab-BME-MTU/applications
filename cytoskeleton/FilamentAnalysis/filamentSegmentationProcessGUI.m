function varargout = filamentSegmentationProcessGUI(varargin)
% filamentSegmentationProcessGUI M-file for filamentSegmentationProcessGUI.fig
%      filamentSegmentationProcessGUI, by itself, creates a new filamentSegmentationProcessGUI or raises the existing
%      singleton*.
%
%      H = filamentSegmentationProcessGUI returns the handle to a new filamentSegmentationProcessGUI or the handle to
%      the existing singleton*. 
%
%      filamentSegmentationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in filamentSegmentationProcessGUI.M with the given input arguments.
%
%      filamentSegmentationProcessGUI('Property','Value',...) creates a new filamentSegmentationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before filamentSegmentationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to filamentSegmentationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filamentSegmentationProcessGUI

% Last Modified by GUIDE v2.5 10-Aug-2014 21:29:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filamentSegmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @filamentSegmentationProcessGUI_OutputFcn, ...
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


% --- Executes just before filamentSegmentationProcessGUI is made visible.
function filamentSegmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',0);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

channelIndex = funParams.ChannelIndex;

% Find any parent process
userData.parentProc = userData.crtPackage.getParent(userData.procID);
if isempty(userData.crtPackage.processes_{userData.procID}) && ~isempty(userData.parentProc)
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,userData.crtPackage.processes_(userData.parentProc)));
    if ~emptyParentProc
        % Intersect channel index with channel index of parent processes
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        for i = userData.parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
    end
   
end

if ~isempty(channelIndex)
    channelString = userData.MD.getChannelPaths(channelIndex);
else
    channelString = {};
end
set(handles.edit_subsample_number,'String',funParams.Sub_Sample_Num);
set(handles.checkbox_outgrowth,'value',funParams.VIF_Outgrowth_Flag);
set(handles.checkbox_rerunwholemovie,'value',funParams.Rerun_WholeMovie);

set(handles.checkbox_nofiguredisruption,'value',funParams.nofiguredisruption);
set(handles.checkbox_savefigures,'value',funParams.savestepfigures);
set(handles.checkbox_showdetailmessage,'value',funParams.showdetailmessages);

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);
set(handles.listbox_selectedChannels,'Value',1);
current_ch_ind = channelIndex(1);

% in the case of single channel or same setting for all channel as in old
% cases
if(length(funParams.StPace_Size)==1)
    ones_array = ones(1,numel(userData.MD.channels_));
    funParams.StPace_Size = funParams.StPace_Size*ones_array;
    funParams.StPatch_Size = funParams.StPatch_Size*ones_array;
    funParams.st_lowerbound_localthresholding = funParams.st_lowerbound_localthresholding*ones_array;
    funParams.IntPace_Size = funParams.IntPace_Size*ones_array;
    funParams.IntPatch_Size = funParams.IntPatch_Size*ones_array;
    funParams.int_lowerbound_localthresholding = funParams.int_lowerbound_localthresholding*ones_array;
    funParams.Cell_Mask_ind = funParams.Cell_Mask_ind*ones_array;
    funParams.Whole_movie_ind = funParams.Whole_movie_ind*ones_array;
    
    Combine_Way = funParams.Combine_Way;
    % if in a cell, unwrap it out as the string
    while(iscell(Combine_Way))
        Combine_Way = Combine_Way{1,1};
    end    
    
    funParams.Combine_Way=cell(1,1);
    
   for iC = 1 : numel(userData.MD.channels_)
        funParams.Combine_Way{iC}= Combine_Way;
    end
    
    funParams.Classifier_Type_ind = funParams.Classifier_Type_ind*ones_array;
    funParams.LengthThreshold = funParams.LengthThreshold*ones_array;
    funParams.CurvatureThreshold = funParams.CurvatureThreshold*ones_array;
    funParams.IternationNumber = funParams.IternationNumber*ones_array;
    funParams.CoefAlpha = funParams.CoefAlpha*ones_array;
    funParams.training_sample_number = funParams.training_sample_number*ones_array;
    funParams.CannyHigherThreshold = funParams.CannyHigherThreshold*ones_array;
    funParams.CannyLowerThreshold = funParams.CannyLowerThreshold*ones_array;
end
    
% with the previous lines of expansion of single value to an array of all
% possible channels, display the first selected channel
set(handles.edit_StPaceSize,'String',funParams.StPace_Size(current_ch_ind));
set(handles.edit_StPatchSize,'String',funParams.StPatch_Size(current_ch_ind));
set(handles.edit_st_lowerbound_localthresholding,'String',funParams.st_lowerbound_localthresholding(current_ch_ind));

set(handles.edit_IntPaceSize,'String',funParams.IntPace_Size(current_ch_ind));
set(handles.edit_IntPatchSize,'String',funParams.IntPatch_Size(current_ch_ind));
set(handles.edit_int_lowerbound_localthresholding,'String',funParams.int_lowerbound_localthresholding(current_ch_ind));

set(handles.popupmenu_cell_mask, 'Value',funParams.Cell_Mask_ind(current_ch_ind));
set(handles.popupmenu_whole_movie,'Value',funParams.Whole_movie_ind(current_ch_ind));


Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
    'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method', 'reserved_for_test'};
for Combine_Way_ind = 1 : 10
    if(strcmp(funParams.Combine_Way{current_ch_ind}, Combine_Way_tag{Combine_Way_ind}))
        set(handles.popupmenu_segmentationbase,'Value',Combine_Way_ind);
        break;
    end
end


% visible or not, show the parameters for the geo based algorithm
% when visible, set the numbers
set(handles.popupmenu_classifier_type,'Value',funParams.Classifier_Type_ind(current_ch_ind));
% current_ch_ind
% funParams.LengthThreshold
set(handles.edit_lengththreshold,'String',funParams.LengthThreshold(current_ch_ind));
set(handles.edit_curvaturethreshold,'String',funParams.CurvatureThreshold(current_ch_ind));
set(handles.edit_IternationNumber,'String',funParams.IternationNumber(current_ch_ind));
set(handles.edit_linear_plane_offset_alpha,'String',funParams.CoefAlpha(current_ch_ind));
set(handles.edit_train_number,'String',funParams.training_sample_number(current_ch_ind));

% show the canny parameters regardless of visibility
set(handles.editCannyHigherThreshold,'String',funParams.CannyHigherThreshold(current_ch_ind));
set(handles.editCannyLowerThreshold,'String',funParams.CannyLowerThreshold(current_ch_ind));


% first set everything as invisible
set(handles.uipanel_threshold_panel,'Visible','off');
set(handles.uipanel_Geo_panel,'Visible','off');
set(handles.uipanel_canny_panel,'Visible','off');



if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_training') ...
        || strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_GM') ...
         || strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_no_GM') ...
      || strcmp(funParams.Combine_Way{current_ch_ind},'reserved_for_test'))
    % with Geo based approaches, use the geo panel and make the
    % thresholding panel invisible
    set(handles.uipanel_threshold_panel,'Visible','off');
    set(handles.uipanel_Geo_panel,'Visible','on');
    set(handles.uipanel_canny_panel,'Visible','off');

    
    if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_training'))
        set(handles.edit_IternationNumber,'Enable','off');
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end
    
    
     
    if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_no_GM'))
        set(handles.edit_IternationNumber,'Enable','off');
        set(handles.edit_IternationNumber,'Value',0');
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end
    
else
    if (strcmp(funParams.Combine_Way{current_ch_ind},'canny_method') )
        set(handles.uipanel_threshold_panel,'Visible','off');
        set(handles.uipanel_Geo_panel,'Visible','off');
        set(handles.uipanel_canny_panel,'Visible','on');

    else
        
        % with thresholding based approaches, use the thresholding panel and make the
        % Geo panel invisible
        set(handles.uipanel_threshold_panel,'Visible','on');
        set(handles.uipanel_Geo_panel,'Visible','off');
        set(handles.uipanel_canny_panel,'Visible','off');
        
        %first enable both sets as default, in the following, details will
        %be set.
        set(handles.edit_StPaceSize,'Enable','on');
        set(handles.edit_StPatchSize,'Enable','on');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    end
end

if (strcmp(funParams.Combine_Way{current_ch_ind},'st_only')...
        ||strcmp(funParams.Combine_Way{current_ch_ind},'st_nms_two')...
        ||strcmp(funParams.Combine_Way{current_ch_ind},'st_nms_only'))
    
    % with st based threhold based approaches, use the threhold panel and make the
    % Geo panel invisible
    set(handles.uipanel_threshold_panel,'Visible','on');
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');
       
    % make the st based parameters enabled, but the intensity based
    % parameters disabled.
    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
end

if (strcmp(funParams.Combine_Way{current_ch_ind},'int_only'))
    % same here, with intensity threhold based approaches, use the threhold panel and make the
    % Geo panel invisible
    
    set(handles.uipanel_threshold_panel,'Visible','on');
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');
       
    
    % make the intensity based parameters enabled, but the st based
    % parameters disabled.
    set(handles.edit_StPaceSize,'Enable','off');
    set(handles.edit_StPatchSize,'Enable','off');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
    set(handles.edit_IntPaceSize,'Enable','on');
    set(handles.edit_IntPatchSize,'Enable','on');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
end

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

processGUI_ApplyFcn_without_close_figure(hObject, eventdata, handles,funParams);


% --- Outputs from this function are returned to the command line.
function varargout = filamentSegmentationProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
currentChannelIndex = channelIndex;

for this_channel_index = currentChannelIndex(:)'
    % if this channel has been specifically signed setting, 
    % skip this channel 
    if(length(funParams.channel_specific)>=this_channel_index)
        if(funParams.channel_specific(this_channel_index)==1)
            continue;
        end
    end
    
    % all possible way of combining
    Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two', ...
        'st_nms_only','geo_based_training','geo_based_GM','geo_based_no_GM','canny_method',...
        'reserved_for_test'};
    
    Combine_Way_ind = get(handles.popupmenu_segmentationbase, 'Value');
    funParams.Combine_Way{this_channel_index}=Combine_Way_tag{Combine_Way_ind};
    
    
    % the patch sizes
    StPace_Size = str2double(get(handles.edit_StPaceSize, 'String'));
    if isnan(StPace_Size) || StPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPace_Size(this_channel_index)=StPace_Size;
    
    StPatch_Size = str2double(get(handles.edit_StPatchSize, 'String'));
    if isnan(StPatch_Size) || StPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPatch_Size(this_channel_index)=StPatch_Size;
    
    
    st_lowerbound_localthresholding = str2double(get(handles.edit_st_lowerbound_localthresholding, 'String'));
    if isnan(st_lowerbound_localthresholding) || st_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.st_lowerbound_localthresholding(this_channel_index)=st_lowerbound_localthresholding;
    
    
    IntPace_Size = str2double(get(handles.edit_IntPaceSize, 'String'));
    if isnan(IntPace_Size) || IntPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPace_Size(this_channel_index)=IntPace_Size;
    
    IntPatch_Size = str2double(get(handles.edit_IntPatchSize, 'String'));
    if isnan(IntPatch_Size) || IntPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPatch_Size(this_channel_index)=IntPatch_Size;
    
    
    int_lowerbound_localthresholding = str2double(get(handles.edit_int_lowerbound_localthresholding, 'String'));
    if isnan(int_lowerbound_localthresholding) || int_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.int_lowerbound_localthresholding(this_channel_index)=int_lowerbound_localthresholding;
    
    Cell_Mask_ind = get(handles.popupmenu_cell_mask, 'Value');
    funParams.Cell_Mask_ind(this_channel_index) = Cell_Mask_ind;
    
    
    Whole_movie_ind = get(handles.popupmenu_whole_movie, 'Value');
    funParams.Whole_movie_ind(this_channel_index) = Whole_movie_ind;
    
    
    if(strcmp(get(handles.uipanel_Geo_panel,'Visible'),'on'))
        Classifier_Type_ind = get(handles.popupmenu_classifier_type,'Value');
        funParams.Classifier_Type_ind(this_channel_index)=Classifier_Type_ind;
        
               
        LengthThreshold = str2double(get(handles.edit_lengththreshold, 'String'));
        if isnan(LengthThreshold) || LengthThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text53,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.LengthThreshold(this_channel_index)=LengthThreshold;
        
        CurvatureThreshold = str2double(get(handles.edit_curvaturethreshold, 'String'));
        if isnan(CurvatureThreshold) || CurvatureThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text52,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.CurvatureThreshold(this_channel_index)=CurvatureThreshold;
        
        CoefAlpha = str2double(get(handles.edit_linear_plane_offset_alpha, 'String'));
        if isnan(CoefAlpha) || CoefAlpha < -1
            errordlg('Please provide a valid input for Alpha in the linear plane classifier','Setting Error','modal');
            return;
        end
        funParams.CoefAlpha(this_channel_index) = CoefAlpha;
        
        training_sample_number = str2double(get(handles.edit_train_number, 'String'));
        if isnan(training_sample_number) || training_sample_number < 0
            errordlg('Please provide a valid input for training sample number','Setting Error','modal');
            return;
        end
        funParams.training_sample_number(this_channel_index) = training_sample_number;
        
        IternationNumber = str2double(get(handles.edit_IternationNumber, 'String'));
        if isnan(IternationNumber) || IternationNumber < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text51,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.IternationNumber(this_channel_index)=IternationNumber;
    end
    
    
    if(strcmp(get(handles.uipanel_canny_panel,'Visible'),'on'))
                       
        CannyHigherThreshold = str2double(get(handles.editCannyHigherThreshold, 'String'));
        if isnan(CannyHigherThreshold) || CannyHigherThreshold < 0  || CannyHigherThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyHigherThreshold(this_channel_index)=CannyHigherThreshold;
        
        CannyLowerThreshold = str2double(get(handles.editCannyLowerThreshold, 'String'));
        if isnan(CannyLowerThreshold) || CannyLowerThreshold < 0  || CannyLowerThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyLowerThreshold(this_channel_index)=CannyLowerThreshold;
    end
    
end

Rerun_WholeMovie = get(handles.checkbox_rerunwholemovie,'value');
funParams.Rerun_WholeMovie = Rerun_WholeMovie;


VIF_Outgrowth_Flag = get(handles.checkbox_outgrowth,'value');
funParams.VIF_Outgrowth_Flag = VIF_Outgrowth_Flag;

nofiguredisruption = get(handles.checkbox_nofiguredisruption,'value');
funParams.nofiguredisruption = nofiguredisruption;

savestepfigures = get(handles.checkbox_savefigures,'value');
funParams.savestepfigures = savestepfigures;

showdetailmessages = get(handles.checkbox_showdetailmessage,'value');
funParams.showdetailmessages = showdetailmessages;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Sub_Sample_Num  = Sub_Sample_Num;
 

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_availableChannels, 'String');

chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
        thresholdValues =zeros(1,numel(chanIndex1));
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
        thresholdValues = [];
end
set(handles.listbox_selectedChannels, 'UserData', chanIndex2);
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues))

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_availableChannels, 'String');
contents2 = get(handles.listbox_selectedChannels, 'String');
id = get(handles.listbox_availableChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_selectedChannels, 'String', contents2, 'Userdata', chanIndex2);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_selectedChannels,'String');
id = get(handles.listbox_selectedChannels,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_selectedChannels, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedChannels,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedChannels,'String',contents);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ishandle(userData.helpFig), delete(userData.helpFig); end
if ishandle(userData.previewFig), delete(userData.previewFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

function edit_StPatchSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StPatchSize as text
%        str2double(get(hObject,'String')) returns contents of edit_StPatchSize as a double


% --- Executes during object creation, after setting all properties.
function edit_StPatchSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_st_lowerbound_localthresholding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_st_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_st_lowerbound_localthresholding as text
%        str2double(get(hObject,'String')) returns contents of edit_st_lowerbound_localthresholding as a double


% --- Executes during object creation, after setting all properties.
function edit_st_lowerbound_localthresholding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_st_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_segmentationbase.
function popupmenu_segmentationbase_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_segmentationbase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_segmentationbase


% The following code is for set the parameter set editing boxes, that corresponds to combine way 
% not choosen as disabled from editing

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
    'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method','reserved_for_test'};

Combine_Way_ind = get(hObject, 'Value');
Combine_Way=Combine_Way_tag{Combine_Way_ind};


set(handles.uipanel_threshold_panel,'Visible','off');
set(handles.uipanel_Geo_panel,'Visible','off');
set(handles.uipanel_canny_panel,'Visible','off');

% see line 122 and on for comments, same thing
if (strcmp(Combine_Way,'st_only')||strcmp(Combine_Way,'st_nms_two')||strcmp(Combine_Way,'st_nms_only'))
     
    set(handles.uipanel_threshold_panel,'Visible','on');    
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');

    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
    
else
    if (strcmp(Combine_Way,'int_only'))
         
    set(handles.uipanel_threshold_panel,'Visible','on');    
    set(handles.uipanel_Geo_panel,'Visible','off');
    
        set(handles.edit_StPaceSize,'Enable','off');
        set(handles.edit_StPatchSize,'Enable','off');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    else
         
    set(handles.uipanel_threshold_panel,'Visible','on');    
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');

        set(handles.edit_StPaceSize,'Enable','on');
        set(handles.edit_StPatchSize,'Enable','on');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    end
end

if (strcmp(Combine_Way,'geo_based_training') || strcmp(Combine_Way,'geo_based_GM')...
        || strcmp(Combine_Way,'geo_based_no_GM')|| strcmp(Combine_Way,'reserved_for_test') )
    
    set(handles.uipanel_threshold_panel,'Visible','off');
    set(handles.uipanel_Geo_panel,'Visible','on');
   set(handles.uipanel_canny_panel,'Visible','off');
    
    if (strcmp(Combine_Way,'geo_based_training'))
        set(handles.edit_IternationNumber,'Enable','off');
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end   
    
      
    if (strcmp(Combine_Way,'geo_based_no_GM'))
        set(handles.edit_IternationNumber,'Enable','off');
        set(handles.edit_IternationNumber,'Value',0');
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end
    
    
    if (strcmp(Combine_Way,'geo_based_GM') ||strcmp(Combine_Way,'geo_based_no_GM') )
          set(handles.popupmenu_classifier_type,'Value',1);  
    end
    
end

if (strcmp(Combine_Way,'canny_method'))
% with canny, no parameter to set for now    
    set(handles.uipanel_threshold_panel,'Visible','off');
    set(handles.uipanel_Geo_panel,'Visible','off');
     set(handles.uipanel_canny_panel,'Visible','on');
  
end

% --- Executes during object creation, after setting all properties.
function popupmenu_segmentationbase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',{'Steerable Filter Results','Intensity','Combine Both St and Int',...
    'St NMS two','St nms only', 'Geo based with training','Geo based with GM',...
    'Geo based without GM','Canny Method',...
    'Reserved for Tests'});

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
    'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method','reserved_for_test'};
% set(hObject,'Value',7);
% for Combine_Way_ind = 1 : 10
%     if(strcmp(funParams.Combine_Way, Combine_Way_tag{Combine_Way_ind}))
%         set(hObject,'Value',Combine_Way_ind);
%         break;
%     end
% end
 

% --- Executes on button press in checkbox_outgrowth.
function checkbox_outgrowth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_outgrowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_outgrowth


% --- Executes on selection change in popupmenu_cell_mask.
function popupmenu_cell_mask_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cell_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cell_mask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cell_mask


% --- Executes during object creation, after setting all properties.
function popupmenu_cell_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cell_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Cell Segmentation in this same channel','Input ROI', ...
    'Cell Segmentation combined from two channels', ...
    'Direct sum of Cell Segmentation channel 1&2','No limitation'});
set(hObject,'Value',4);



function edit_subsample_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subsample_number as text
%        str2double(get(hObject,'String')) returns contents of edit_subsample_number as a double


% --- Executes during object creation, after setting all properties.
function edit_subsample_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IntPaceSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IntPaceSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IntPaceSize as text
%        str2double(get(hObject,'String')) returns contents of edit_IntPaceSize as a double


% --- Executes during object creation, after setting all properties.
function edit_IntPaceSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IntPaceSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IntPatchSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IntPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IntPatchSize as text
%        str2double(get(hObject,'String')) returns contents of edit_IntPatchSize as a double


% --- Executes during object creation, after setting all properties.
function edit_IntPatchSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IntPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_int_lowerbound_localthresholding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_int_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_int_lowerbound_localthresholding as text
%        str2double(get(hObject,'String')) returns contents of edit_int_lowerbound_localthresholding as a double


% --- Executes during object creation, after setting all properties.
function edit_int_lowerbound_localthresholding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_int_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_traing_nms.
function pushbutton_traing_nms_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_traing_nms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
%  
% if isempty(userData.MD)    
%     userData_main.MD = MD;    
% else    
%     userData_main.MD = cat(2, userData_main.MD, MD);    
% end

funParams.F_classifier = nms_classifier_train_moviedata(userData.MD, channelIndex);

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in pushbutton_load_nms_classifier.
function pushbutton_load_nms_classifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_nms_classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This code direct user to copy the trained .mat file, and check if it is
% currently in place.
% Yes or not, it is assumed the user have done so or will do so, hence set
% the file name for classifier as if it is in place.

userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
F_classifer_train_output = cell(1,max(channelIndex));

nPackage = length(userData.MD.packages_);

indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(userData.MD.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

nProcesses = length(userData.MD.processes_);

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(userData.MD.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

FilamentSegmentationProcessOutputDir  = [userData.MD.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
if (~exist(FilamentSegmentationProcessOutputDir,'dir'))
    mkdir(FilamentSegmentationProcessOutputDir);
end

for iChannel = channelIndex
    
%     Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  [FilamentSegmentationProcessOutputDir, '/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
       
    F_classifer_train_output{iChannel} = [FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'];
    
    if(~exist([FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'],'file'))
       f = showinfowindow(['Please copy the trained classifier as F_classifer_channel.mat in ...\FilamentAnalysisPackage\FilamentSegmentation\Channel',...
        num2str(iChannel),' for the trained movie, and paste it at corresponding position for current movie.','  Currently cannot find the copied F\_classifer\_channel.mat for current movie, if you have not done so, please copy and double check.'], 'Warning');
    else
       f = showinfowindow(['Please copy the trained classifier as F_classifer_channel.mat in ...\FilamentAnalysisPackage\FilamentSegmentation\Channel',...
        num2str(iChannel),' for the trained movie, and paste it at corresponding position for current movie.','  Currently it is detected that the .mat file is copied.'], 'Message');
        
    end    
end

funParams.F_classifier = F_classifer_train_output;

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



% --- Executes during object creation, after setting all properties.
function pushbutton_traing_nms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_traing_nms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_deletetraining.
function pushbutton_deletetraining_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletetraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
 
funParams.F_classifier = cell(1,max(channelIndex));

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_showdetailmessage.
function checkbox_showdetailmessage_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showdetailmessage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showdetailmessage


% --- Executes on button press in checkbox_nofiguredisruption.
function checkbox_nofiguredisruption_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nofiguredisruption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nofiguredisruption


% --- Executes on button press in checkbox_savefigures.
function checkbox_savefigures_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savefigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_savefigures



function edit_IternationNumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IternationNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IternationNumber as text
%        str2double(get(hObject,'String')) returns contents of edit_IternationNumber as a double


% --- Executes during object creation, after setting all properties.
function edit_IternationNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IternationNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_classifier_type.
function popupmenu_classifier_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_classifier_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_classifier_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_classifier_type

Classifier_Type_tag = {'linear','SVM'};

Classifier_Type_ind = get(hObject, 'Value');
Classifier_Type=Classifier_Type_tag{Classifier_Type_ind};


% --- Executes during object creation, after setting all properties.
function popupmenu_classifier_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_classifier_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'linear', 'SVM'});

Combine_Way_tag = {'linear','SVM'};
set(hObject,'Value',1);


function edit_curvaturethreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_curvaturethreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_curvaturethreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_curvaturethreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_curvaturethreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_curvaturethreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lengththreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lengththreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lengththreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_lengththreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_lengththreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lengththreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_linear_plane_offset_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_linear_plane_offset_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_linear_plane_offset_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_linear_plane_offset_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_linear_plane_offset_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_linear_plane_offset_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_train_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_train_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_train_number as text
%        str2double(get(hObject,'String')) returns contents of edit_train_number as a double


% --- Executes during object creation, after setting all properties.
function edit_train_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_train_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_whole_movie.
function popupmenu_whole_movie_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_whole_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_whole_movie contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_whole_movie


% --- Executes during object creation, after setting all properties.
function popupmenu_whole_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_whole_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Completely', 'Half-half','None'});

Whole_movie_flag = {'Completely', 'Half-half','None'};
set(hObject,'Value',2);


% --- Executes on button press in pushbutton_para_this_channel.
function pushbutton_para_this_channel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_para_this_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply to this channel' button
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

channelIndex = get(handles.listbox_selectedChannels, 'Userdata');

funParams.ChannelIndex = channelIndex;

blued_channel = get(handles.listbox_selectedChannels, 'Value');
currentChannelIndex = channelIndex(blued_channel(1));

for this_channel_index = currentChannelIndex(:)
    % Mark 1 if this channel has been specifically signed setting
    
    funParams.channel_specific(this_channel_index)=1;
   
    Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
        'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method','reserved_for_test'};
    Combine_Way_ind = get(handles.popupmenu_segmentationbase, 'Value');
    funParams.Combine_Way{this_channel_index}=Combine_Way_tag{Combine_Way_ind};
    
    StPace_Size = str2double(get(handles.edit_StPaceSize, 'String'));
    if isnan(StPace_Size) || StPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPace_Size(this_channel_index)=StPace_Size;
    
    StPatch_Size = str2double(get(handles.edit_StPatchSize, 'String'));
    if isnan(StPatch_Size) || StPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPatch_Size(this_channel_index)=StPatch_Size;
    
    
    st_lowerbound_localthresholding = str2double(get(handles.edit_st_lowerbound_localthresholding, 'String'));
    if isnan(st_lowerbound_localthresholding) || st_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.st_lowerbound_localthresholding(this_channel_index)=st_lowerbound_localthresholding;
    
    
    IntPace_Size = str2double(get(handles.edit_IntPaceSize, 'String'));
    if isnan(IntPace_Size) || IntPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPace_Size(this_channel_index)=IntPace_Size;
    
    IntPatch_Size = str2double(get(handles.edit_IntPatchSize, 'String'));
    if isnan(IntPatch_Size) || IntPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPatch_Size(this_channel_index)=IntPatch_Size;
    
    
    int_lowerbound_localthresholding = str2double(get(handles.edit_int_lowerbound_localthresholding, 'String'));
    if isnan(int_lowerbound_localthresholding) || int_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.int_lowerbound_localthresholding(this_channel_index)=int_lowerbound_localthresholding;
    
    Cell_Mask_ind = get(handles.popupmenu_cell_mask, 'Value');
    funParams.Cell_Mask_ind(this_channel_index) = Cell_Mask_ind;
    
    
    Whole_movie_ind = get(handles.popupmenu_whole_movie, 'Value');
    funParams.Whole_movie_ind(this_channel_index) = Whole_movie_ind;
    
    
    if(strcmp(get(handles.uipanel_Geo_panel,'Visible'),'on'))
        Classifier_Type_ind = get(handles.popupmenu_classifier_type,'Value');
        funParams.Classifier_Type_ind(this_channel_index)=Classifier_Type_ind;
        
        LengthThreshold = str2double(get(handles.edit_lengththreshold, 'String'));
        if isnan(LengthThreshold) || LengthThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text53,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.LengthThreshold(this_channel_index)=LengthThreshold;
        
        CurvatureThreshold = str2double(get(handles.edit_curvaturethreshold, 'String'));
        if isnan(CurvatureThreshold) || CurvatureThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text52,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.CurvatureThreshold(this_channel_index)=CurvatureThreshold;
        
        CoefAlpha = str2double(get(handles.edit_linear_plane_offset_alpha, 'String'));
        if isnan(CoefAlpha) || CoefAlpha < -1
            errordlg('Please provide a valid input for Alpha in the linear plane classifier','Setting Error','modal');
            return;
        end
        funParams.CoefAlpha(this_channel_index) = CoefAlpha;
        
        training_sample_number = str2double(get(handles.edit_train_number, 'String'));
        if isnan(training_sample_number) || training_sample_number < 0
            errordlg('Please provide a valid input for training sample number','Setting Error','modal');
            return;
        end
        funParams.training_sample_number(this_channel_index) = training_sample_number;
        
        IternationNumber = str2double(get(handles.edit_IternationNumber, 'String'));
        if isnan(IternationNumber) || IternationNumber < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text51,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.IternationNumber(this_channel_index)=IternationNumber;
    end
    
    
    if(strcmp(get(handles.uipanel_canny_panel,'Visible'),'on'))
                       
        CannyHigherThreshold = str2double(get(handles.editCannyHigherThreshold, 'String'));
        if isnan(CannyHigherThreshold) || CannyHigherThreshold < 0  || CannyHigherThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyHigherThreshold(this_channel_index)=CannyHigherThreshold;
        
        CannyLowerThreshold = str2double(get(handles.editCannyLowerThreshold, 'String'));
        if isnan(CannyLowerThreshold) || CannyLowerThreshold < 0  || CannyLowerThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyLowerThreshold(this_channel_index)=CannyLowerThreshold;
    end
    
    msgbox(['Setting assigned to channel', num2str(this_channel_index)]);
end

Rerun_WholeMovie = get(handles.checkbox_rerunwholemovie,'value');
funParams.Rerun_WholeMovie = Rerun_WholeMovie;

VIF_Outgrowth_Flag = get(handles.checkbox_outgrowth,'value');
funParams.VIF_Outgrowth_Flag = VIF_Outgrowth_Flag;

nofiguredisruption = get(handles.checkbox_nofiguredisruption,'value');
funParams.nofiguredisruption = nofiguredisruption;

savestepfigures = get(handles.checkbox_savefigures,'value');
funParams.savestepfigures = savestepfigures;

showdetailmessages = get(handles.checkbox_showdetailmessage,'value');
funParams.showdetailmessages = showdetailmessages;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Sub_Sample_Num  = Sub_Sample_Num;
    
processGUI_ApplyFcn_without_close_figure(hObject, eventdata, handles,funParams);

% --- Executes on button press in pushbutton_para_all_selected_channels.
function pushbutton_para_all_selected_channels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_para_all_selected_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply to all selected channel' button
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
currentChannelIndex = channelIndex;
funParams.LengthThreshold=[];
  
for this_channel_index = currentChannelIndex(:)'
    % even if this channel has been specifically signed setting, 
    % re-assign acorrding to the current one.
    funParams.channel_specific(this_channel_index)=0;
    Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
        'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method','reserved_for_test'};
    Combine_Way_ind = get(handles.popupmenu_segmentationbase, 'Value');
    funParams.Combine_Way{this_channel_index} = Combine_Way_tag{Combine_Way_ind};
    
    StPace_Size = str2double(get(handles.edit_StPaceSize, 'String'));
    if isnan(StPace_Size) || StPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPace_Size(this_channel_index)=StPace_Size;
    
    StPatch_Size = str2double(get(handles.edit_StPatchSize, 'String'));
    if isnan(StPatch_Size) || StPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.StPatch_Size(this_channel_index)=StPatch_Size;
    
    
    st_lowerbound_localthresholding = str2double(get(handles.edit_st_lowerbound_localthresholding, 'String'));
    if isnan(st_lowerbound_localthresholding) || st_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.st_lowerbound_localthresholding(this_channel_index)=st_lowerbound_localthresholding;
    
    
    IntPace_Size = str2double(get(handles.edit_IntPaceSize, 'String'));
    if isnan(IntPace_Size) || IntPace_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPace_Size(this_channel_index)=IntPace_Size;
    
    IntPatch_Size = str2double(get(handles.edit_IntPatchSize, 'String'));
    if isnan(IntPatch_Size) || IntPatch_Size < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.IntPatch_Size(this_channel_index)=IntPatch_Size;
    
    
    int_lowerbound_localthresholding = str2double(get(handles.edit_int_lowerbound_localthresholding, 'String'));
    if isnan(int_lowerbound_localthresholding) || int_lowerbound_localthresholding < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
        return;
    end
    funParams.int_lowerbound_localthresholding(this_channel_index)=int_lowerbound_localthresholding;
    
    Cell_Mask_ind = get(handles.popupmenu_cell_mask, 'Value');
    funParams.Cell_Mask_ind(this_channel_index) = Cell_Mask_ind;
    
    
    Whole_movie_ind = get(handles.popupmenu_whole_movie, 'Value');
    funParams.Whole_movie_ind(this_channel_index) = Whole_movie_ind;
    
   if(strcmp(get(handles.uipanel_Geo_panel,'Visible'),'on'))
        Classifier_Type_ind = get(handles.popupmenu_classifier_type,'Value');
        funParams.Classifier_Type_ind(this_channel_index)=Classifier_Type_ind;
        
        LengthThreshold = str2double(get(handles.edit_lengththreshold, 'String'));
        if isnan(LengthThreshold) || LengthThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text53,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.LengthThreshold(this_channel_index)=LengthThreshold;
        
        CurvatureThreshold = str2double(get(handles.edit_curvaturethreshold, 'String'));
        if isnan(CurvatureThreshold) || CurvatureThreshold < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text52,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.CurvatureThreshold(this_channel_index)=CurvatureThreshold;
        
        CoefAlpha = str2double(get(handles.edit_linear_plane_offset_alpha, 'String'));
        if isnan(CoefAlpha) || CoefAlpha < -1
            errordlg('Please provide a valid input for Alpha in the linear plane classifier','Setting Error','modal');
            return;
        end
        funParams.CoefAlpha(this_channel_index) = CoefAlpha;
        
        training_sample_number = str2double(get(handles.edit_train_number, 'String'));
        if isnan(training_sample_number) || training_sample_number < 0
            errordlg('Please provide a valid input for training sample number','Setting Error','modal');
            return;
        end
        funParams.training_sample_number(this_channel_index) = training_sample_number;
        
        IternationNumber = str2double(get(handles.edit_IternationNumber, 'String'));
        if isnan(IternationNumber) || IternationNumber < 0
            errordlg(['Please provide a valid input for '''...
                get(handles.text51,'String') '''.'],'Setting Error','modal');
            return;
        end
        funParams.IternationNumber(this_channel_index)=IternationNumber;
   end
    
   
   if(strcmp(get(handles.uipanel_canny_panel,'Visible'),'on'))
                       
        CannyHigherThreshold = str2double(get(handles.editCannyHigherThreshold, 'String'));
        if isnan(CannyHigherThreshold) || CannyHigherThreshold < 0  || CannyHigherThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyHigherThreshold(this_channel_index)=CannyHigherThreshold;
        
        CannyLowerThreshold = str2double(get(handles.editCannyLowerThreshold, 'String'));
        if isnan(CannyLowerThreshold) || CannyLowerThreshold < 0  || CannyLowerThreshold >100
            errordlg('Please provide a valid input for Canny Higher Threshold',...
                'Setting Error','modal');
            return;
        end
        funParams.CannyLowerThreshold(this_channel_index)=CannyLowerThreshold;
    end
end

Rerun_WholeMovie = get(handles.checkbox_rerunwholemovie,'value');
funParams.Rerun_WholeMovie = Rerun_WholeMovie;


VIF_Outgrowth_Flag = get(handles.checkbox_outgrowth,'value');
funParams.VIF_Outgrowth_Flag = VIF_Outgrowth_Flag;

nofiguredisruption = get(handles.checkbox_nofiguredisruption,'value');
funParams.nofiguredisruption = nofiguredisruption;

savestepfigures = get(handles.checkbox_savefigures,'value');
funParams.savestepfigures = savestepfigures;

showdetailmessages = get(handles.checkbox_showdetailmessage,'value');
funParams.showdetailmessages = showdetailmessages;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Sub_Sample_Num  = Sub_Sample_Num;
 
msgbox('Setting assigned to all selected channels');

processGUI_ApplyFcn_without_close_figure(hObject, eventdata, handles,funParams);

% 
% % --- Executes on key press with focus on listbox_selectedChannels and none of its controls.
% function listbox_selectedChannels_KeyPressFcn(hObject, eventdata, handles)
% % hObject    handle to listbox_selectedChannels (see GCBO)
% % eventdata  structure with the following fields (see UICONTROL)
% %	Key: name of the key that was pressed, in lower case
% %	Character: character interpretation of the key(s) that was pressed
% %	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% 
% 
% % --- If Enable == 'on', executes on mouse press in 5 pixel border.
% % --- Otherwise, executes on mouse press in 5 pixel border or over listbox_selectedChannels.
% function listbox_selectedChannels_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to listbox_selectedChannels (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% --- Executes on selection change in listbox_selectedChannels.
function listbox_selectedChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectedChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_selectedChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selectedChannels

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

channelIndex = funParams.ChannelIndex;
blued_channel = get(handles.listbox_selectedChannels, 'Value');
currentChannelIndex = channelIndex(blued_channel(1));

set(handles.edit_subsample_number,'String',funParams.Sub_Sample_Num);
set(handles.checkbox_outgrowth,'value',funParams.VIF_Outgrowth_Flag);
set(handles.checkbox_rerunwholemovie,'value',funParams.Rerun_WholeMovie);

set(handles.checkbox_nofiguredisruption,'value',funParams.nofiguredisruption);
set(handles.checkbox_savefigures,'value',funParams.savestepfigures);
set(handles.checkbox_showdetailmessage,'value',funParams.showdetailmessages);

current_ch_ind = currentChannelIndex(1);

%% comment these since in the opening this has been run, so this code below
% should never be active.

% % in the case of single channel or same setting for all channel as in old
% % cases
% if(length(funParams.StPace_Size)==1)
%     ones_array = ones(1,numel(userData.MD.channels_));
%     funParams.StPace_Size = funParams.StPace_Size*ones_array;
%     funParams.StPatch_Size = funParams.StPatch_Size*ones_array;
%     funParams.st_lowerbound_localthresholding = funParams.st_lowerbound_localthresholding*ones_array;
%     funParams.IntPace_Size = funParams.IntPace_Size*ones_array;
%     funParams.IntPatch_Size = funParams.IntPatch_Size*ones_array;
%     funParams.int_lowerbound_localthresholding = funParams.int_lowerbound_localthresholding*ones_array;
%     funParams.Cell_Mask_ind = funParams.Cell_Mask_ind*ones_array;
%     funParams.Whole_movie_ind = funParams.Whole_movie_ind*ones_array;
%     
%     Combine_Way = funParams.Combine_Way;
%     for iC = 1 : numel(userData.MD.channels_)
%         funParams.Combine_Way{iC}= Combine_Way;
%     end
%     
%     funParams.Classifier_Type_ind = funParams.Classifier_Type_ind*ones_array;
%     funParams.LengthThreshold = funParams.LengthThreshold*ones_array;
%     funParams.CurvatureThreshold = funParams.CurvatureThreshold*ones_array;
%     funParams.IternationNumber = funParams.IternationNumber*ones_array;
%     funParams.CoefAlpha = funParams.CoefAlpha*ones_array;
%     funParams.training_sample_number = funParams.training_sample_number*ones_array;
% end
%     
% with the previous lines of expansion of single value to an array of all
% possible channels, display the first selected channel
set(handles.edit_StPaceSize,'String',funParams.StPace_Size(current_ch_ind));
set(handles.edit_StPatchSize,'String',funParams.StPatch_Size(current_ch_ind));
set(handles.edit_st_lowerbound_localthresholding,'String',funParams.st_lowerbound_localthresholding(current_ch_ind));

set(handles.edit_IntPaceSize,'String',funParams.IntPace_Size(current_ch_ind));
set(handles.edit_IntPatchSize,'String',funParams.IntPatch_Size(current_ch_ind));
set(handles.edit_int_lowerbound_localthresholding,'String',funParams.int_lowerbound_localthresholding(current_ch_ind));

set(handles.popupmenu_cell_mask, 'Value',funParams.Cell_Mask_ind(current_ch_ind));
set(handles.popupmenu_whole_movie,'Value',funParams.Whole_movie_ind(current_ch_ind));

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only',...
    'geo_based_training','geo_based_GM','geo_based_no_GM','canny_method','reserved_for_test'};
for Combine_Way_ind = 1 : 10
    if(strcmp(funParams.Combine_Way{current_ch_ind}, Combine_Way_tag{Combine_Way_ind}))
        set(handles.popupmenu_segmentationbase,'Value',Combine_Way_ind);
        break;
    end
end

% visible or not, show the parameters for the geo based algorithm
% when visible, set the numbers
set(handles.popupmenu_classifier_type,'Value',funParams.Classifier_Type_ind(current_ch_ind));
set(handles.edit_lengththreshold,'String',funParams.LengthThreshold(current_ch_ind));
set(handles.edit_curvaturethreshold,'String',funParams.CurvatureThreshold(current_ch_ind));
set(handles.edit_IternationNumber,'String',funParams.IternationNumber(current_ch_ind));
set(handles.edit_linear_plane_offset_alpha,'String',funParams.CoefAlpha(current_ch_ind));
set(handles.edit_train_number,'String',funParams.training_sample_number(current_ch_ind));

% first set everything as invisible
set(handles.uipanel_threshold_panel,'Visible','off');
set(handles.uipanel_Geo_panel,'Visible','off');
set(handles.uipanel_canny_panel,'Visible','off');

if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_training') ...
        || strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_GM') ...
         || strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_no_GM') ...
        || strcmp(funParams.Combine_Way{current_ch_ind},'reserved_for_test'))
    % with Geo based approaches, use the geo panel and make the
    % thresholding panel invisible
    set(handles.uipanel_threshold_panel,'Visible','off');
    set(handles.uipanel_Geo_panel,'Visible','on');
     set(handles.uipanel_canny_panel,'Visible','off');
   
    if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_training'))
        set(handles.edit_IternationNumber,'Enable','off');
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end
    
    if (strcmp(funParams.Combine_Way{current_ch_ind},'geo_based_no_GM'))
        set(handles.edit_IternationNumber,'Enable','off');
         set(handles.edit_IternationNumber,'Value',0);
    else
        set(handles.edit_IternationNumber,'Enable','on');
    end
    
else
    
    if (strcmp(funParams.Combine_Way{current_ch_ind},'canny_method') )
        set(handles.uipanel_threshold_panel,'Visible','off');
        set(handles.uipanel_Geo_panel,'Visible','off');
        set(handles.uipanel_canny_panel,'Visible','on');

    else
    
    % with thresholding based approaches, use the thresholding panel and make the
    % Geo panel invisible
    set(handles.uipanel_threshold_panel,'Visible','on');
    set(handles.uipanel_Geo_panel,'Visible','off');
     set(handles.uipanel_canny_panel,'Visible','off');

    %first enable both sets as default, in the following, details will
    %be set.
    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','on');
    set(handles.edit_IntPatchSize,'Enable','on');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    end
end

if (strcmp(funParams.Combine_Way{current_ch_ind},'st_only')||strcmp(funParams.Combine_Way{current_ch_ind},'st_nms_two')||strcmp(funParams.Combine_Way{current_ch_ind},'st_nms_only'))
    
    % with st based threhold based approaches, use the threhold panel and make the
    % Geo panel invisible
    set(handles.uipanel_threshold_panel,'Visible','on');
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');

    % make the st based parameters enabled, but the intensity based
    % parameters disabled.
    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
end

if (strcmp(funParams.Combine_Way{current_ch_ind},'int_only'))
    % same here, with intensity threhold based approaches, use the threhold panel and make the
    % Geo panel invisible
    
    set(handles.uipanel_threshold_panel,'Visible','on');
    set(handles.uipanel_Geo_panel,'Visible','off');
    set(handles.uipanel_canny_panel,'Visible','off');

    % make the intensity based parameters enabled, but the st based
    % parameters disabled.
    set(handles.edit_StPaceSize,'Enable','off');
    set(handles.edit_StPatchSize,'Enable','off');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
    set(handles.edit_IntPaceSize,'Enable','on');
    set(handles.edit_IntPatchSize,'Enable','on');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
end


% % since this is for display not input, no need for update
%  processGUI_ApplyFcn_without_close_figure(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_rerunwholemovie.
function checkbox_rerunwholemovie_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rerunwholemovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rerunwholemovie



function editCannyLowerThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editCannyLowerThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCannyLowerThreshold as text
%        str2double(get(hObject,'String')) returns contents of editCannyLowerThreshold as a double


% --- Executes during object creation, after setting all properties.
function editCannyLowerThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCannyLowerThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCannyHigherThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editCannyHigherThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCannyHigherThreshold as text
%        str2double(get(hObject,'String')) returns contents of editCannyHigherThreshold as a double


% --- Executes during object creation, after setting all properties.
function editCannyHigherThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCannyHigherThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
