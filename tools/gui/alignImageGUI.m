function varargout = alignImageGUI(varargin)
% ALIGNIMAGEGUI M-file for alignImageGUI.fig
%      ALIGNIMAGEGUI, by itself, creates a new ALIGNIMAGEGUI or raises the existing
%      singleton*.
%
%      H = ALIGNIMAGEGUI returns the handle to a new ALIGNIMAGEGUI or the handle to
%      the existing singleton*.
%
%      ALIGNIMAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNIMAGEGUI.M with the given input arguments.
%
%      ALIGNIMAGEGUI('Property','Value',...) creates a new ALIGNIMAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignImageGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignImageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Author: Lin Ji, May 01, 2006.

% Edit the above text to modify the response to help alignImageGUI

% Last Modified by GUIDE v2.5 17-Apr-2006 17:10:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alignImageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @alignImageGUI_OutputFcn, ...
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


% --- Executes just before alignImageGUI is made visible.
function alignImageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignImageGUI (see VARARGIN)

% Choose default command line output for alignImageGUI
handles.output = hObject;

handles = setDefPar(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alignImageGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alignImageGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function handles = setDefPar(handles)

handles.refImgFile        = '';
handles.alignImgPathList  = {'---- Empty ----'};
handles.numAlignImgPaths  = 0;
handles.totalNumImages    = 0;
handles.imgWidth          = 0;
handles.imgHeight         = 0;
handles.alignFrame        = 0;
handles.modelImgPathID    = 0;
handles.firstAlignImgFile = {};
handles.markerROI         = [];
handles.maxShift          = 50;
handles.alignShift        = [];

handles.alignImgFigH = [];
handles.refImgFigH   = [];

set(handles.refImageEdit,'String',handles.refImgFile);
set(handles.alignFrameEdit,'String',handles.alignFrame);
set(handles.alignImgPathMenu,'String',handles.alignImgPathList);
set(handles.maxShiftEdit,'String',handles.maxShift);
set(handles.firstImgEdit,'String','');

% --- Executes on button press in refImagePB.
function refImagePB_Callback(hObject, eventdata, handles)
% hObject    handle to refImagePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[refImgFile,pathname,filterindex] = uigetfile({'*.tif';'*.jpg';'*.png';'*.gif'}, ...
   'Select reference image');

if filterindex == 0
   %It is cancelled by the user.
   return;
end

handles.refImgFile = [pathname filesep refImgFile];
set(handles.refImageEdit,'String',handles.refImgFile);

guidata(hObject,handles);


function refImageEdit_Callback(hObject, eventdata, handles)
% hObject    handle to refImageEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refImageEdit as text
%        str2double(get(hObject,'String')) returns contents of refImageEdit as a double

refImgFile = get(handles.refImageEdit,'String');
if exist(refImgFile,'file') == 2
   handles.refImgFile = refImgFile;
   guidata(hObject,handles);
else
   set(handles.refImageEdit,'String',handles.refImgFile);
end


% --- Executes during object creation, after setting all properties.
function refImageEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refImageEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxShiftEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refImageEdit as text
%        str2double(get(hObject,'String')) returns contents of refImageEdit as a double

maxShiftStr = get(handles.maxShiftEdit,'String');

handles.maxShift = str2num(maxShiftStr);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function maxShiftEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function alignFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to alignFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refImageEdit as text
%        str2double(get(hObject,'String')) returns contents of refImageEdit as a double

alignFrameStr = get(handles.alignFrameEdit,'String');

handles.alignFrame = str2num(alignFrameStr);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function alignFrameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in displayPB.
function displayPB_Callback(hObject, eventdata, handles)
% hObject    handle to displayPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dispPath = handles.modelImgPathID;
if dispPath == 0
   dispPath = 1;
end
alignImgFileList = getFileStackNames([handles.alignImgPathList{dispPath} filesep ...
   handles.firstAlignImgFile{dispPath}]);

if handles.alignFrame == 0
   dispImgNo = 1;
else
   dispImgNo = handles.alignFrame;
end

dispImg = imread(alignImgFileList{dispImgNo});

if isempty(handles.alignImgFigH) || ~ishandle(handles.alignImgFigH)
   handles.alignImgFigH = figure;
end

figure(handles.alignImgFigH); hold off;
imshow(dispImg,[]); hold on;
if ~isempty(handles.markerROI)
   plot(handles.markerROI(:,1),handles.markerROI(:,2),'y-.');
end

if dispImgNo <= size(handles.alignShift,1)
   alignShift = handles.alignShift(dispImgNo,:);
elseif isempty(handles.alignShift)
   alignShift = handles.alignShift;
else
   alignShift = handles.alignShift(1,:);
end

if ~isempty(alignShift) 
   if ~isnan(alignShift(1))
      plot(handles.markerROI(:,1)+alignShift(1),handles.markerROI(:,2)+alignShift(2),'w');
      title(sprintf('Aligning Shift (pixels): (%d,%d)',alignShift(1),alignShift(2)));
   else
      title(sprintf('Not Aligned'));
   end
else
   title(sprintf('Not Aligned'));
end

guidata(hObject,handles);

% --- Executes on button press in markerROIPB.
function markerROIPB_Callback(hObject, eventdata, handles)
% hObject    handle to markerROIPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.refImgFile)
   return;
end

refImg = imread(handles.refImgFile);

if isempty(handles.refImgFigH) || ~ishandle(handles.refImgFigH)
   handles.refImgFigH = figure;
end
figure(handles.refImgFigH); hold off;
imshow(refImg,[]); hold on;
if ~isempty(handles.markerROI)
   plot(handles.markerROI(:,1),handles.markerROI(:,2),'w-.');
end

[bw,xi,yi] = roipoly;
handles.markerROI = [xi yi];
plot(handles.markerROI(:,1),handles.markerROI(:,2),'w');

%Indicate that the ROI has been selected.
set(handles.markerROIPB,'ForegroundColor','green');

guidata(hObject,handles);

% --- Executes on button press in startAlignPB.
function startAlignPB_Callback(hObject, eventdata, handles)
% hObject    handle to startAlignPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

alignFrame = handles.alignFrame;

if alignFrame == 0
   handles.alignShift = alignImage(handles.refImgFile,handles.alignImgPathList, ...
      handles.firstAlignImgFile,'modelChannel',handles.modelImgPathID, ...
      'alignFrame',handles.alignFrame,'markerROI',handles.markerROI, ...
      'maxXShift',handles.maxShift,'maxYShift',handles.maxShift);
else
   handles.alignShift(alignFrame,:) = alignImage(handles.refImgFile,handles.alignImgPathList, ...
      handles.firstAlignImgFile,'modelChannel',handles.modelImgPathID, ...
      'alignFrame',handles.alignFrame,'markerROI',handles.markerROI, ...
      'maxXShift',handles.maxShift,'maxYShift',handles.maxShift);
end

guidata(hObject,handles);

if alignFrame == 0
   %For monitoring process status.
   L = length(num2str(handles.totalNumImages));
   strForm = sprintf('%%.%dd',L);

   %Align all images and save them.
   for k = 1:length(handles.alignImgPathList)
      %Create the directory for saving the aligned images.
      alignImgPath = handles.alignImgPathList{k};

      %Remove 'filesep' from the end.
      while strcmp(alignImgPath(end),filesep)
         alignImgPath(end) = '';
      end

      saveAlignImgPath = [alignImgPath '_align'];
      if ~isdir(saveAlignImgPath)
         success = mkdir(saveAlignImgPath);
         if ~success
            fprintf(1,'Trouble making directory. Aligned images not saved.\n');
            return;
         end
      end

      alignImgFileList = getFileStackNames([handles.alignImgPathList{k} filesep ...
         handles.firstAlignImgFile{k}]);

      %Align all images according to the calculated 'alignShift' and save it.
      imgWidth  = handles.imgWidth;
      imgHeight = handles.imgHeight;

      fprintf(1,'Save aligned image: ');
      for ii = 1:length(alignImgFileList)
         procStr = sprintf(strForm,ii);
         fprintf(1,procStr);

         xOffSet = handles.alignShift(ii,1);
         yOffSet = handles.alignShift(ii,2);

         img2Align  = imread(alignImgFileList{ii});
         alignedImg = img2Align;
         alignedImg(:) = 0;
         if xOffSet > 0
            alignedImg(:,1:end-xOffSet) = img2Align(:,1+xOffSet:end);
         else
            alignedImg(:,1-xOffSet:end) = img2Align(:,1:end+xOffSet);
         end
         if yOffSet > 0
            alignedImg(1:end-yOffSet,:) = alignedImg(1+yOffSet:end,:);
         else
            alignedImg(1-yOffSet:end,:) = alignedImg(1:end+yOffSet,:);
         end

         %Save the aligned image.
         [path,body,no,ext] = getFilenameBody(alignImgFileList{ii});
         saveAlignImgFile = [saveAlignImgPath filesep body '_align' no ext];
         imwrite(alignedImg,saveAlignImgFile,'tiff');

         for jj = 1:length(procStr)
            fprintf(1,'\b');
         end
      end
      fprintf(1,[strForm '\n'],handles.totalNumImages);
   end
end

% --- Executes during object creation, after setting all properties.
function startAlignPB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startAlignPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in alignImgPathMenu.
function alignImgPathMenu_Callback(hObject, eventdata, handles)
% hObject    handle to alignImgPathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns alignImgPathMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alignImgPathMenu

if handles.numAlignImgPaths == 0
   return;
else
   handles.modelImgPathID = get(handles.alignImgPathMenu,'Value');
   set(handles.firstImgEdit,'String',handles.firstAlignImgFile{handles.modelImgPathID});
end

% --- Executes during object creation, after setting all properties.
function alignImgPathMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignImgPathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function firstImgEdit_Callback(hObject, eventdata, handles)
% hObject    handle to firstImgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstImgEdit as text
%        str2double(get(hObject,'String')) returns contents of firstImgEdit as a double


% --- Executes during object creation, after setting all properties.
function firstImgEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstImgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsePB.
function browsePB_Callback(hObject, eventdata, handles)
% hObject    handle to browsePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[firstImgFile,pathName,filterindex] = uigetfile({'*.tif';'*.jpg';'*.png';'*.gif'}, ...
   'Add one image channel for aligning by selecting the first image');

if filterindex == 0
   %Job is cancelled.
   return;
end

%First check whether the selected image path exist in handles.alignImgPathList.
k = 1;
while k <= length(handles.numAlignImgPaths)
   if samdir(pathName,handles.alignImgPathList{k})
      return;
   else
      k = k+1;
   end
end

alignImgFileList = getFileStackNames([pathName filesep firstImgFile]);
if handles.totalNumImages == 0
   handles.totalNumImages = length(alignImgFileList);
else
   handles.totalNumImages = min(handles.totalNumImages,length(alignImgFileList));
end
handles.alignShift = NaN*ones(handles.totalNumImages,2);

firstImg = imread(alignImgFileList{1});
handles.imgHeight = size(firstImg,1);
handles.imgWidth  = size(firstImg,2);

if handles.numAlignImgPaths == 0
   handles.alignImgPathList{1} = pathName;
else
   handles.alignImgPathList{end+1} = pathName;
end
handles.firstAlignImgFile{end+1} = firstImgFile;
handles.numAlignImgPaths = handles.numAlignImgPaths+1;

%Make the added image path the model path.
handles.modelImgPathID = handles.numAlignImgPaths;

set(handles.alignImgPathMenu,'String',handles.alignImgPathList);
set(handles.alignImgPathMenu,'Value',handles.modelImgPathID);
set(handles.firstImgEdit,'String',handles.firstAlignImgFile{handles.modelImgPathID});

guidata(hObject,handles);


% --- Executes on button press in deletePB.
function deletePB_Callback(hObject, eventdata, handles)
% hObject    handle to deletePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.numAlignImgPaths == 0
   return;
end

handles.alignImgPathList(handles.modelImgPathID)  = [];
handles.firstAlignImgFile(handles.modelImgPathID) = [];

if handles.modelImgPathID == handles.numAlignImgPaths
   handles.modelImgPathID = handles.modelImgPathID-1;
end
handles.numAlignImgPaths = handles.numAlignImgPaths-1;

if handles.numAlignImgPaths == 0
   handles.alignImgPathList = {'---- Empty ----'};
   set(handles.alignImgPathMenu,'Value',1);
   set(handles.firstImgEdit,'String','');
else
   set(handles.alignImgPathMenu,'Value',handles.modelImgPathID);
   set(handles.firstImgEdit,'String',handles.firstAlignImgFile{handles.modelImgPathID});
end

set(handles.alignImgPathMenu,'String',handles.alignImgPathList);

guidata(hObject,handles);

