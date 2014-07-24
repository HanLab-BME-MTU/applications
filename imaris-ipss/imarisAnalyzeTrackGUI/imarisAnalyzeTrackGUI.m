function varargout = imarisAnalyzeTrackGUI(varargin)
% IATGUI M-file for imarisAnalyzeTrackGUI.fig
%      IATGUI, by itself, creates a new IATGUI or raises the existing
%      singleton*.
%
%      H = IATGUI returns the handle to a new IATGUI or the handle to
%      the existing singleton*.
%
%      IATGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IATGUI.M with the given input arguments.
%
%      IATGUI('Property','Value',...) creates a new IATGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imarisAnalyzeTrackGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imarisAnalyzeTrackGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help imarisAnalyzeTrackGUI

% Last Modified by GUIDE v2.5 29-Nov-2004 15:00:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imarisAnalyzeTrackGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @imarisAnalyzeTrackGUI_OutputFcn, ...
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


% --- Executes just before imarisAnalyzeTrackGUI is made visible.
function imarisAnalyzeTrackGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imarisAnalyzeTrackGUI (see VARARGIN)

% Choose default command line output for imarisAnalyzeTrackGUI
handles.output = hObject;

% empty varargin - incorrect call
if isempty(varargin)
    error('imarisAnalyzeTrackGUI needs tracks as input arguments. Run imarisAnalyzeTrack.')
end

% get data and store
handles.tracks = varargin{1};
% count tracks so that we know what options to offer
numTracks = length(handles.tracks);

%collect handles
handles.selectorHandles = [handles.IAT_select1_PD];

handles.select1Handles = [handles.IAT_track11_txt,...
        handles.IAT_track11_PD;...
        handles.IAT_track12_txt,...
        handles.IAT_track12_PD];


%get data for pulldown menus
PD_data = iatGUI_PD_data;

% remove all entries for which we do not have enough tracks
requiredNumOfTracks = cat(1,PD_data{:,2});
PD_data(requiredNumOfTracks>numTracks,:) = [];
% store
handles.PD_data = PD_data;

% set track names
trackNames = cell(numTracks+1,1);
trackNames{1} = 'please select...';
[trackNames{2:end}] = deal(handles.tracks.name);
set(handles.select1Handles(:,2),'Value',1,'String',trackNames);


%set labelMenus and hide tag-PD (none -> no )
set(handles.selectorHandles,'String',handles.PD_data(:,1));
set([handles.select1Handles],'Visible','off');

%init fields
handles.plotFigureH = [];

%remember last selection of data PD
handles.lastSelected = 1;


% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = imarisAnalyzeTrackGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in IAT_select1_PD.
function IAT_select1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to IAT_select1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%activate necessary tag-PD's
selectedItem = get(hObject,'Value');

necessaryTracks=handles.PD_data{selectedItem, 2};

if necessaryTracks == -1
    %user selected a separator
    set(hObject,'Value',handles.lastSelected(1));
else
    set(handles.select1Handles(:),'Visible','off');
    set(handles.select1Handles(1:necessaryTracks,:),'Visible','on');
    handles.lastSelected(1) = selectedItem;
    
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function IAT_select1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IAT_select1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in IAT_track11_PD.
function IAT_track11_PD_Callback(hObject, eventdata, handles)
% hObject    handle to IAT_track11_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function IAT_track11_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IAT_track11_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in IAT_track12_PD.
function IAT_track12_PD_Callback(hObject, eventdata, handles)
% hObject    handle to IAT_track12_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns IAT_track12_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IAT_track12_PD


% --- Executes during object creation, after setting all properties.
function IAT_track12_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IAT_track12_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in IAT_close_PB.
function IAT_close_PB_Callback(hObject, eventdata, handles)
% hObject    handle to IAT_close_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ans = questdlg('Do you really want to close this panel and all figures created with it?',...
    '','Yes','No','Yes');
switch ans
    case 'Yes'
        delete(handles.plotFigureH);
        delete(handles.iatGUI);
    case 'No'
        % do nothing
        return
end
        

% --- Executes on button press in IAT_plot_PB.
function IAT_plot_PB_Callback(hObject, eventdata, handles)
% hObject    handle to IAT_plot_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 1. check that all selections have been made
% 2. choose plotting function
% 3. read and pass all the necessary data

plotSelection = get(handles.IAT_select1_PD,'Value');
trackSelection = get(handles.select1Handles(:,2),'Value');

numRequiredTracks = handles.PD_data{plotSelection,2};

% checks
if numRequiredTracks < 1
    h= errordlg('please select something to plot');
    uiwait(h)
    return
end

for i=1:numRequiredTracks
    if trackSelection{i} == 1
        h = errordlg('please select all necessary tags');
        uiwait(h);
        return
    end
end

% choose what to plot
plotFunction = handles.PD_data{plotSelection,4};

switch plotFunction
%     1: distance to origin
%     2: distance between tracks
%     3: displacement (along track)
%     4: displacement (relative to track)
%     5: MSQD relative to origin
%     6: MSQD relative to track
%     7: MSQD relative to track, 1D
    
    case 2 % distance to origin
        
        % we need: pos1, fill in zeros for pos2
        trackNum1 = trackSelection{1}-1;
        trackName1 = handles.tracks(trackNum1).name;
        trackName2 = 'origin';
        spots1 = handles.tracks(trackNum1).spotXYZ;
        spots2 = zeros(size(spots1));
        xData = handles.tracks(trackNum1).spotT;
        
        % calculate distance
        yData = iatGUI_distance(spots1, spots2);
        
        % legend text
        legendTxt = ['Distance ' trackName1 ' - ' trackName2];
        labelX = 'Time';
        labelY = 'Distance';
        
        % sigmas
        doSigma = 0;
        
    case 3 % distance between tracks
        
        % we need: pos1, pos2
        trackNum1 = trackSelection{1}-1;
        trackNum2 = trackSelection{2}-1;
        trackName1 = handles.tracks(trackNum1).name;
        trackName2 = handles.tracks(trackNum2).name;
        spots1 = handles.tracks(trackNum1).spotXYZ;
        spots2 = handles.tracks(trackNum2).spotXYZ;
        
        % make sure we have only time where there are two spots
        time1 = handles.tracks(trackNum1).spotT;
        time2 = handles.tracks(trackNum2).spotT;
        
        [xData, idx1, idx2] = intersect(time1, time2);
        spots1 = spots1(idx1,:);
        spots2 = spots2(idx2,:);
        
        % calculate distance
        yData = iatGUI_distance(spots1, spots2);
        
        % legend text
        legendTxt = ['Distance ' trackName1 ' - ' trackName2];    
        labelX = 'Time';
        labelY = 'Distance';
        
        % sigmas
        doSigma = 0;
        
    case 4 % displacement (along track)
        
         % we need: pos1, zeros for pos 2 then delta
         trackNum1 = trackSelection{1}-1;
        trackName1 = handles.tracks(trackNum1).name;
        spots1 = handles.tracks(trackNum1).spotXYZ;
        spots2 = zeros(size(spots1));
        
        xData = handles.tracks(trackNum1).spotT;
        xData = xData(1:end-1);
        
        % calculate distance
        yData = iatGUI_displacement(spots1, spots2);
        
        % legend text
        legendTxt = ['Displacement of ' trackName1];    
        labelX = 'Time';
        labelY = 'Displacement';
        
        % sigmas
        doSigma = 0;
        
        
    case 5 % displacement (relative to track)
        
        % we need: pos1, pos2, then diff&delta
        trackNum1 = trackSelection{1}-1;
        trackNum2 = trackSelection{2}-1;
        trackName1 = handles.tracks(trackNum1).name;
        trackName2 = handles.tracks(trackNum2).name;
        spots1 = handles.tracks(trackNum1).spotXYZ;
        spots2 = handles.tracks(trackNum2).spotXYZ;
        
        % make sure we have only time where there are two spots
        time1 = handles.tracks(trackNum1).spotT;
        time2 = handles.tracks(trackNum2).spotT;
        
        [xData, idx1, idx2] = intersect(time1, time2);
        xData = xData(1:end-1);
        spots1 = spots1(idx1,:);
        spots2 = spots2(idx2,:);
        
        % calculate distance
        yData = iatGUI_displacement(spots1, spots2);
        
        % legend text
        legendTxt = ['Displacement of ' trackName1 ' relative to ' trackName2];    
        labelX = 'Time';
        labelY = 'Displacement';
        
        % sigmas
        doSigma = 0;
        
        
    case 6 % MSQD relative to origin
        
        % need: pos1,  then calc msqd
        
        trackNum1 = trackSelection{1}-1;
        trackName1 = handles.tracks(trackNum1).name;
        spTmp1 = handles.tracks(trackNum1).spotXYZ;
        
        % make sure we have only time where there are two spots
        time1 = handles.tracks(trackNum1).spotT;
        
        % we have to fill spots with NaNs where there is nothing to keep a
        % constant time interval. Time index will always be with increment
        % 1 (or so I hope)
        numTime = max(time1);
        allTime = 1:numTime;
        
        spots1 = repmat(NaN,numTime,3);
        % allIdx tell where the goodTimes are placed in the allTime matrix
        [dummy, dummy, allIdx] = intersect(allTime, time1);
        
        % take goodTime - spots and place into allTime matrices
        spots1(allIdx,:) = spTmp1(time1,:);
        
        % no known covariances, therefore use zeros
        covMat = zeros(3,3,numTime);
        
        % assing input for msqd
        positions(1).coordinates = spots1;
        positions(1).covariances = covMat;
        doOneD = 0;
        
        % calculate msqd and assign
        msqd = iatGUI_MSQD(positions,doOneD);
        yData = msqd(:,1);
        ySigma = msqd(:,2);
        xData = allTime(1:length(yData)); % get correct time!
        
        % legend text
        legendTxt = ['MSQD of ' trackName1 ' relative to origin'];    
        labelX = 'Timelag';
        labelY = 'MSQD';
        
        % sigma
        doSigma = 1;
        
    case 7 % MSQD relative to track
        
        % need: pos1, pos2, then calc msqd
        
        trackNum1 = trackSelection{1}-1;
        trackNum2 = trackSelection{2}-1;
        trackName1 = handles.tracks(trackNum1).name;
        trackName2 = handles.tracks(trackNum2).name;
        spTmp1 = handles.tracks(trackNum1).spotXYZ;
        spTmp2 = handles.tracks(trackNum2).spotXYZ;
        
        % make sure we have only time where there are two spots
        time1 = handles.tracks(trackNum1).spotT;
        time2 = handles.tracks(trackNum2).spotT;
        
        % we have to fill spots with NaNs where there is nothing to keep a
        % constant time interval. Time index will always be with increment
        % 1 (or so I hope)
        numTime = max(max(time1),max(time2));
        allTime = 1:numTime;
        
        [spots1,spots2] = deal(repmat(NaN,numTime,3));
        [goodTime, idx1, idx2] = intersect(time1, time2);
        % allIdx tell where the goodTimes are placed in the allTime matrix
        [dummy, dummy, allIdx] = intersect(allTime, goodTime);
        
        % take goodTime - spots and place into allTime matrices
        spots1(allIdx,:) = spTmp1(idx1,:);
        spots2(allIdx,:) = spTmp2(idx2,:);
        
        % no known covariances, therefore use zeros
        covMat = zeros(3,3,numTime);
        
        % assing input for msqd
        positions(1).coordinates = spots1;
        positions(1).covariances = covMat;
        positions(2).coordinates = spots2;
        positions(2).covariances = covMat;
        doOneD = 0;
        
        % calculate msqd and assign
        msqd = iatGUI_MSQD(positions,doOneD);
        yData = msqd(:,1);
        ySigma = msqd(:,2);
        xData = allTime(1:length(yData)); % get correct time!
        
        % legend text
        legendTxt = ['MSQD of ' trackName1 ' relative to ' trackName2];    
        labelX = 'Timelag';
        labelY = 'MSQD';
        
        % sigma
        doSigma = 1;
        
        
    case 8 % MSQD relative to track, 1D
        
        % need: pos1, pos2, then calc msqd, doOneD = 1
        
        trackNum1 = trackSelection{1}-1;
        trackNum2 = trackSelection{2}-1;
        trackName1 = handles.tracks(trackNum1).name;
        trackName2 = handles.tracks(trackNum2).name;
        spTmp1 = handles.tracks(trackNum1).spotXYZ;
        spTmp2 = handles.tracks(trackNum2).spotXYZ;
        
        % make sure we have only time where there are two spots
        time1 = handles.tracks(trackNum1).spotT;
        time2 = handles.tracks(trackNum2).spotT;
        
        % we have to fill spots with NaNs where there is nothing to keep a
        % constant time interval. Time index will always be with increment
        % 1 (or so I hope)
        numTime = max(max(time1),max(time2));
        allTime = 1:numTime;
        
        [spots1,spots2] = deal(repmat(NaN,numTime,3));
        [goodTime, idx1, idx2] = intersect(time1, time2);
        % allIdx tell where the goodTimes are placed in the allTime matrix
        [dummy, dummy, allIdx] = intersect(allTime, goodTime);
        
        % take goodTime - spots and place into allTime matrices
        spots1(allIdx,:) = spTmp1(idx1,:);
        spots2(allIdx,:) = spTmp2(idx2,:);
        
        % no known covariances, therefore use zeros
        covMat = zeros(3,3,numTime);
        
        % assing input for msqd
        positions(1).coordinates = spots1;
        positions(1).covariances = covMat;
        positions(2).coordinates = spots2;
        positions(2).covariances = covMat;
        doOneD = 1;
        
        % calculate msqd and assign
        msqd = iatGUI_MSQD(positions,doOneD);
        yData = msqd(:,1);
        ySigma = msqd(:,2);
        xData = allTime(1:length(yData)); % get correct time!
        
        % legend text
        legendTxt = ['1D MSQD of ' trackName1 ' relative to ' trackName2];    
        labelX = 'Timelag';
        labelY = 'MSQD';
        
        % sigma
        doSigma = 1;
        
    otherwise
        
        h = errordlg('sorry, function is not implemented yet');
        uiwait(h);
        return
        
end


% plot, write legend
fh = figure;
ph = plot(xData,yData,'d-');
xlabel(labelX);
ylabel(labelY);
if doSigma
    myErrorbar(xData,yData,ySigma);
end
legend(ph, legendTxt);

% store figure handle
handles.plotFigureH = [handles.plotFigureH; fh];
guidata(hObject,handles);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SUBFUNCTIONS

function distance = iatGUI_distance(spots1, spots2)
distance = normList(spots1-spots2);

function displacement = iatGUI_displacement(spots1, spots2)
pos = spots1 - spots2;
displacement = normList(pos(2:end,:)-pos(1:end-1,:));

function mSqDisp = iatGUI_MSQD(positions,doOneD)
mSqDisp = meanSquaredDisplacement(positions, [], doOneD);
        
        
    

