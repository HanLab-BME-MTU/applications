function varargout = guiSim(varargin)
% GUISIM MATLAB code for guiSim.fig
%      GUISIM, by itself, creates a new GUISIM or raises the existing
%      singleton*.
%
%      H = GUISIM returns the handle to a new GUISIM or the handle to
%      the existing singleton*.
%
%      GUISIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUISIM.M with the given input arguments.
%
%      GUISIM('Property','Value',...) creates a new GUISIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiSim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiSim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiSim

% Last Modified by GUIDE v2.5 25-Jan-2012 14:23:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiSim_OpeningFcn, ...
                   'gui_OutputFcn',  @guiSim_OutputFcn, ...
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


% --- Executes just before guiSim is made visible.
function guiSim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiSim (see VARARGIN)

% Choose default command line output for guiSim
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiSim wait for user response (see UIRESUME)
% uiwait(handles.figure1);

data = guidata(gcf);

data.size = [1000,1000,250];
data.stdDev = [5,5,5];
data.samplingType = 'regular';
data.samplingDensity = 1/10;
data.newCurve = [];
data.models = {};
data.out = Data();
data.pdf = [];
data.stdDevCorrFactor = 1;
data.nNoisePoints = 0;
data.roiCenter = [0 0 0];

data.controlPointSize = 7;
data.mousePosXY = [0,0];
data.mousePosXZ = [0,0];
data.isDraggingXY = [];
data.isDraggingXZ = [];
data.maxCureDeg = 50;

global controlKey;
controlKey = false;

% Setup GUI
set(handles.figure1,'UserData',false); 
set(handles.figure1,'Units','Pixels');
set(handles.figure1,'WindowButtonMotionFcn',{@mouseCursorFcn,handles});
set(handles.figure1,'KeyPressFcn',{@keyPressCallback,handles});
set(handles.figure1,'KeyReleaseFcn',{@keyReleaseCallback,handles})
set(handles.axesXY,'Units','Pixels');
set(handles.axesXZ,'Units','Pixels');
set(handles.axesXY,'Box','on');
set(handles.axesXZ,'Box','on');

% Display instructions
instructions(1) = {'Ctrl: Keep pressed and click to add new curves'};
instructions(2) = {'Enter: Abort adding control points'};
instructions(3) = {'Shift: Move control points'};
instructions(4) = {'D: Delete curves'};
text(-400,0,instructions);

% Update GUI
data.initialAxesSize = get(handles.axesXY,'Position');
guidata(gcf,data);

updateAxes(handles,data.initialAxesSize);
updateGUIfields(handles);


function updateAxes(handles,initialAxesSize)
data = guidata(gcf);
scaleFactorXY = min(initialAxesSize(3:4)./data.size(1:2));

previewSizeXY = scaleFactorXY * data.size(1:2);
set(handles.axesXY,'Position',[initialAxesSize(1:2),previewSizeXY]);
axis(handles.axesXY,[data.roiCenter(1)-data.size(1)/2 data.roiCenter(1)+data.size(1)/2 data.roiCenter(2)-data.size(2)/2 data.roiCenter(2)+data.size(2)/2]);

posAxesXZ = get(handles.axesXZ,'Position');

previewSizeXZ = scaleFactorXY * data.size([1,3]);
set(handles.axesXZ,'Position',[posAxesXZ(1:2),previewSizeXZ]);
axis(handles.axesXZ,[data.roiCenter(1)-data.size(1)/2 data.roiCenter(1)+data.size(1)/2 data.roiCenter(3)-data.size(3)/2 data.roiCenter(3)+data.size(3)/2]);


function updateGUIfields(handles)
data = guidata(gcf);
set(handles.edit12,'String',num2str(data.stdDev(1)));
set(handles.edit13,'String',num2str(data.stdDev(2)));
set(handles.edit14,'String',num2str(data.stdDev(3)));

set(handles.edit15,'String',num2str(data.samplingDensity));

set(handles.edit6,'String',num2str(data.size(1)));
set(handles.edit7,'String',num2str(data.size(2)));
set(handles.edit8,'String',num2str(data.size(3)));

set(handles.edit16,'String',num2str(data.stdDevCorrFactor));
set(handles.edit17,'String',num2str(data.nNoisePoints));

if ~isempty(data.pdf)
    set(handles.radiobutton7,'Visible','on');
else
    set(handles.radiobutton7,'Visible','off');
end
set(handles.radiobutton4,'Value',1);


function keyPressCallback(src,evnt,handles)
data = guidata(gcf);

global controlKey;

k = evnt.Key;

% Add control point to the new curve
if strcmp(k,'control')
    if ~controlKey
        controlKey = true;
        appendPoint2NewCurve(handles);
        addNewCurve(handles);
        controlKey = false;
    end
end
% Drag point
if strcmp(k,'shift')
    nModels = numel(data.models);
    minDistXY = [];
    minDistXZ = [];

    for m=1:nModels
        cP = data.models{m};
        for c=1:size(cP,1)
            distXY = norm(cP(c,1:2)-data.mousePosXY);
            distXZ = norm(cP(c,[1,3])-data.mousePosXZ);

            if distXY <= data.controlPointSize && (isempty(minDistXY) || distXY < minDistXY)
                minDistXY = distXY;
                data.isDraggingXY = [m,c];
                % disp('hit XY')
            elseif distXZ <= data.controlPointSize && (isempty(minDistXZ) || distXZ < minDistXZ)
                minDistXZ = distXZ;
                data.isDraggingXZ = [m,c];
                % disp('hit XZ')
            end
        end
    end
    guidata(gcf,data);
end
% Delete model
if strcmp(k,'d')
    nModels = numel(data.models);
    minDistXY = [];
    minDistXZ = [];
    model2delete = [];
    
    for m=1:nModels
        cP = data.models{m};
        for c=1:size(cP,1)
            distXY = norm(cP(c,1:2)-data.mousePosXY);
            distXZ = norm(cP(c,[1,3])-data.mousePosXZ);

            if distXY <= data.controlPointSize && (isempty(minDistXY) || distXY < minDistXY)
                minDistXY = distXY;
                model2delete = m;
                % disp('hit XY')
            elseif distXZ <= data.controlPointSize && (isempty(minDistXZ) || distXZ < minDistXZ)
                minDistXZ = distXZ;
                model2delete = m;
                % disp('hit XZ')
            end
        end
    end
    % Delete model
    if ~isempty(model2delete)
        idx = true(numel(data.models),1);
        idx(model2delete) = false;
        data.models = data.models(idx);
        guidata(gcf,data);
        drawModels(handles);
    end
end


function keyReleaseCallback(src,evnt,handles)
data = guidata(gcf);

global controlKey;
k = evnt.Key;

if strcmp(k,'control')
    controlKey = false;
    if ispc
        import java.awt.Robot;
        import java.awt.event.*;
        keyboard = Robot;
        keyboard.keyPress(KeyEvent.VK_ENTER);
    end
end
if strcmp(k,'shift')
    data.isDraggingXY = [];
    data.isDraggingXZ = [];
    guidata(gcf,data);
end


function appendPoint2NewCurve(handles)
data = guidata(gcf);
global controlKey;

while controlKey
    [x,y,button] = ginput(1);
    if ~controlKey
        break;
    end
    if ~isempty(button)
        data.newCurve(end+1,1:3) = [x,y,data.size(3)/2];
        guidata(gcf,data);
        drawModels(handles);
    end
end


function addNewCurve(handles)
data = guidata(gcf);
% Check if array contains at least two control points
if size(data.newCurve,1) > 1
    data.models(end+1) = {data.newCurve(1:min(end,data.maxCureDeg),:)};
end
data.newCurve = [];
guidata(gcf,data);
drawModels(handles);


function drawModels(handles)
data = guidata(gcf);

nModels = numel(data.models);
cla(handles.axesXY);
cla(handles.axesXZ);

n = 20;

for m=1:nModels
    cP = data.models{m};
    drawCurve(cP,n,handles,data);
end
% Plot new curve
cP = data.newCurve(1:min(size(data.newCurve,1),data.maxCureDeg),:);
if size(cP,1) > 1
    drawCurve(cP,n,handles,data);
end


function drawCurve(cP,n,handles,data)
% Plot the Bezier curve
t = linspace(0,1,n)';
XYZ = renderBezier(cP,t);
pXY = plot(handles.axesXY,XYZ(:,1),XYZ(:,2));
set(pXY,'Color','red','LineWidth',2)
hold(handles.axesXY,'on');
pXZ = plot(handles.axesXZ,XYZ(:,1),XYZ(:,3));
set(pXZ,'Color','red','LineWidth',2)
hold(handles.axesXZ,'on');
axis(handles.axesXY,[0 data.size(1) 0 data.size(2)]);
axis(handles.axesXZ,[0 data.size(1) 0 data.size(3)]);
% Plot the control points
for c=1:size(cP,1)
    plot(handles.axesXY,cP(c,1),cP(c,2),'ko','LineWidth',2,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerSize',data.controlPointSize);
    plot(handles.axesXZ,cP(c,1),cP(c,3),'ko','LineWidth',2,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerSize',data.controlPointSize);
end
% Plot the tangents
if size(cP,1) > 2
    plot(handles.axesXY,cP([1,2],1),cP([1,2],2));
    plot(handles.axesXY,cP([end,end-1],1),cP([end,end-1],2));
    plot(handles.axesXZ,cP([1,2],1),cP([1,2],3));
    plot(handles.axesXZ,cP([end,end-1],1),cP([end,end-1],3));
end


function mouseCursorFcn(hObject,eventdata,handles)
data = guidata(gcf);
% for XY
newPointXY = get(hObject,'CurrentPoint');
posXY = get(handles.axesXY,'Position');
inAxesXY = (newPointXY(1,1)>=posXY(1)) && ...
    (newPointXY(1,1)<posXY(1)+posXY(3)) && ...
    (newPointXY(1,2)>posXY(2)) && ...
    (newPointXY(1,2)<posXY(2)+posXY(4));
% for XZ
newPointXZ = get(hObject,'CurrentPoint');
posXZ = get(handles.axesXZ,'Position');
inAxesXZ = (newPointXZ(1,1)>=posXZ(1)) && ...
    (newPointXZ(1,1)<posXZ(1)+posXZ(3)) && ...
    (newPointXZ(1,2)>posXZ(2)) && ...
    (newPointXZ(1,2)<posXZ(2)+posXZ(4));
if inAxesXY
    posXY = get(handles.axesXY,'CurrentPoint');
    data.mousePosXY = posXY(1,1:2);
    guidata(gcf,data);
    if ~isempty(data.isDraggingXY)
        data.models{data.isDraggingXY(1)}(data.isDraggingXY(2),1:2) = data.mousePosXY;
        guidata(gcf,data);
        drawModels(handles);
    end
elseif inAxesXZ
    posXZ = get(handles.axesXZ,'CurrentPoint');
    data.mousePosXZ = posXZ(1,1:2);
    guidata(gcf,data);
    if ~isempty(data.isDraggingXZ)
        data.models{data.isDraggingXZ(1)}(data.isDraggingXZ(2),[1,3]) = data.mousePosXZ;
        guidata(gcf,data);
        drawModels(handles);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = guiSim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
disp('e6')
disp('Width ...')
data = guidata(gcf);
string = get(hObject,'String');
value = str2double(string);
if ~isnan(value)
    data.size(1) = value;
end
guidata(gcf,data);
updateAxes(handles,data.initialAxesSize);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
disp('e7')
disp('Height ...')
data = guidata(gcf);
string = get(hObject,'String');
value = str2double(string);
if ~isnan(value)
    data.size(2) = value;
end
guidata(gcf,data);
updateAxes(handles,data.initialAxesSize);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
disp('e8')
disp('Depth ...')
data = guidata(gcf);
string = get(hObject,'String');
value = str2double(string);
if ~isnan(value)
    data.size(3) = value;
end
guidata(gcf,data);
updateAxes(handles,data.initialAxesSize);


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Save ...')
data = guidata(gcf);
% Init data
data.out = Data();
% Get file name
[file,path] = uiputfile('*.d.dat','Save Data As');
if path ~= 0
    % Generate points
    sim = Simulation(data.out);
    switch data.samplingType
        case 'regular'
            sim.setSamplingToRegular();
            sim.modelToPoints(data.models,data.stdDev*data.stdDevCorrFactor,data.samplingDensity);
        case 'random'
            sim.setSamplingToRandom();
            sim.modelToPoints(data.models,data.stdDev*data.stdDevCorrFactor,data.samplingDensity);
        case 'pdf'
            sim.modelToPointsFromPDF(data.models,data.stdDev*data.stdDevCorrFactor,data.pdf);
    end

    sim.setDomain(data.roiCenter-data.size/2,data.size);
    sim.addRandomNoise(data.nNoisePoints);
    data.nNoisePoints
    if data.stdDevCorrFactor ~= 1
        warndlg(sprintf('The generated data has a standard deviation that was corrected with stdDevCorrFactor = %.2f to match the original data set. However, the corresponding data.error array was not modified.',data.stdDevCorrFactor),'Notification!')
    end
    % Write data object
    data.out.roiSize = [data.size(1:2) 0];
    if data.stdDev(3) ~= 0
        data.out.error = repmat(data.stdDev,size(data.out.points,1),1);
    else
        data.out.error = repmat([data.stdDev(1:2) 1],size(data.out.points,1),1);
        disp('Main: 2D data: Error in Z is set to 1!');
    end
    data.out.simModelBezCP = data.models;
    % mkdir(path,file(1:end-6))
    % data.out.save([path file(1:end-6) '\' file]);
    data.out.save([path file]);
    guidata(gcf,data);
end


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
disp('e12')
disp('x ...')
data = guidata(gcf);
string = get(handles.edit12,'String');
value = str2double(string);
if ~isnan(value)
    data.stdDev(1) = value;
end
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
disp('13')
disp('y ...')
data = guidata(gcf);
string = get(handles.edit13,'String');
value = str2double(string);
if ~isnan(value)
    data.stdDev(2) = value;
end
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
disp('e14')
disp('z ...')
data = guidata(gcf);
string = get(handles.edit14,'String');
value = str2double(string);
if ~isnan(value)
    data.stdDev(3) = value;
end
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Load ...')
data = guidata(gcf);
% Open file
[file,path] = uigetfile('*.dat','Open Data ...');
if path ~= 0
    % Read data object
    in = Data.load([path file]);
    data.out = Data();
    % Warning if no simulation data available
    if ~isempty(in.simModelBezCP) || ~isempty(in.modelBezCP)
        % Make model choice
        if isempty(in.simModelBezCP)
            choice = questdlg('', ...
                'Select a model ...', ...
                'Reconstruction','No, thank you!','Reconstruction');
        elseif isempty(in.modelBezCP)
            choice = questdlg('', ...
                'Select a model ...', ...
                'Ground Truth','No, thank you!','Ground Truth');
        else
            choice = questdlg('', ...
                'Select a model ...', ...
                'Reconstruction','Ground Truth','No, thank you!','Ground Truth');
        end
        % Handle response
        data.stdDevCorrFactor = 1;
        switch choice
            case 'Reconstruction'
                mod = in.modelBezCP;
                ana = Analysis(in);
                [~,data.pdf] = ana.spacings('off');
                data.stdDevCorrFactor = ana.width('off')/sqrt(2);
                data.nNoisePoints = numel(in.nullCluster);
                delete(ana);
                disp('Probability density function loaded!');
            case 'Ground Truth'
                mod = in.simModelBezCP;
                data.pdf = [];
                data.stdDevCorrFactor = 1;
                data.nNoisePoints = 0;
            case 'No, thank you!'
                mod = [];
        end
        if ~isempty(mod)
            % Copy fields
            data.size = in.roiSize;
            data.roiCenter = [0 0 0];
            if in.roiSize(3) == 0
                data.size(3) = 1.25/0.9*diff(quantile(in.points(:,3),[.05 .95])); % Add 10 percent
                data.roiCenter(3) = sum(quantile(in.points(:,3),[.05 .95]))/2;
            end
            data.stdDev = in.error(1,:);
            data.samplingType = 'regular';
            data.models = mod;
            length = cellfun(@lengthBezier,mod);
            totalModelLength = sum(length);
            data.samplingDensity = 1/round(totalModelLength/in.nPoints);
            data.newCurve = [];
            guidata(gcf,data);
            drawModels(handles);
            updateGUIfields(handles);
            updateAxes(handles,data.initialAxesSize);
        else
            disp('WARNING: No model data will be loaded!');
        end
    else
        disp('WARNING: No model data found!');
    end
end


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
disp('rb13')
data = guidata(gcf);
if get(handles.radiobutton4,'Value')
    data.samplingType = 'regular';
    set(handles.edit15,'Visible','on');
    set(handles.text12,'Visible','on');
    disp('regular')
elseif get(handles.radiobutton5,'Value')
    data.samplingType = 'random';
    set(handles.edit15,'Visible','on');
    set(handles.text12,'Visible','on');
    disp('random')
elseif get(handles.radiobutton7,'Value')
    data.samplingType = 'pdf';
    set(handles.edit15,'Visible','off');
    set(handles.text12,'Visible','off');
    disp('pdf')
end
guidata(gcf,data);


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
disp('e15')
disp('density ...')
data = guidata(gcf);
string = get(handles.edit15,'String');
value = str2double(string);
if ~isnan(value)
    data.samplingDensity = value;
end
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
disp('e16')
data = guidata(gcf);
data.stdDevCorrFactor = str2double(get(hObject,'String'));
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
disp('e17')
data = guidata(gcf);
data.nNoisePoints = str2double(get(hObject,'String'));
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
