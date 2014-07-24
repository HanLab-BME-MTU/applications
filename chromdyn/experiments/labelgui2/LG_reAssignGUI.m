function varargout = LG_reAssignGUI(varargin)
% LG_REASSIGNGUI M-file for LG_reAssignGUI.fig
%      LG_REASSIGNGUI, by itself, creates a new LG_REASSIGNGUI or raises the existing
%      singleton*.
%
%      H = LG_REASSIGNGUI returns the handle to a new LG_REASSIGNGUI or the handle to
%      the existing singleton*.
%
%      LG_REASSIGNGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LG_REASSIGNGUI.M with the given input arguments.
%
%      LG_REASSIGNGUI('Property','Value',...) creates a new LG_REASSIGNGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LG_reAssignGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LG_reAssignGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help LG_reAssignGUI

% Last Modified by GUIDE v2.5 18-Nov-2005 18:38:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LG_reAssignGUI_OpeningFcn, ...
    'gui_OutputFcn',  @LG_reAssignGUI_OutputFcn, ...
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


% --- Executes just before LG_reAssignGUI is made visible.
function LG_reAssignGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LG_reAssignGUI (see VARARGIN)

% assign output - need to do here in case we quit
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% INIT GUI

% collect handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;
reAssignH = hObject;
reAssignHandles = get(reAssignH);

% get current time
currentTime = LG_getCurrentTime;
idlist = movieWindowHandles.idlist;
linklist = idlist(currentTime).linklist;


% if no linklist for current time: Don't quit. Instead, open the GUI
% without tags.
if isempty(linklist)
    noTags = true;
else
    noTags = false;
end

% gui-dependent defaults
def_textPos = [2,0,12,1.7];
def_pd1Pos = [16,0,15,1.7];
def_pd2Pos = [33,0,15,1.7];
def_minY = 12;

if noTags
    nGoodSpots = 0;
else
% read idlist
idlistData = movieWindowHandles.idlistData;

% find good spots/tags - all but estimated and single occurences. Correct for
% fusions
goodSpotsIdx = ...
    find(linklist(:,5) ~= 3 & linklist(:,3) ~= 3 & linklist(:,3) ~= 1);
nGoodSpots = length(goodSpotsIdx);

% good tags include fusion tags!
goodTagIdx = find(linklist(:,5) ~= 3);
end


% in case there weren't any good spots, we could delete the frame here
% (load full idlist!), and quit again.
% [idlist,status] = LG_deleteFrame(idlist,currentTime,goodTimes,recalc);

% set height of reAssignGUI, and adjust position.
if ~isempty(naviHandles.positions.LG_reAssignGUI)
    oldPosition = naviHandles.positions.LG_reAssignGUI;
else
    oldPosition = reAssignHandles.Position;
end
% newHeight is the minimum (12) plus one line for every tag plus one title
% line plus  a char for the border
newHeight = def_minY + 2 * (nGoodSpots + 1) + 1;
oldHeight = oldPosition(4);
% take into account the shift in y position
newPosition = oldPosition;
newPosition([2,4]) = [oldPosition(2) + oldHeight - newHeight, newHeight];
% assign new position
set(reAssignH,'Position',newPosition)


if ~noTags

% prepare writing the lines of the GUI
lc = idlistData.labelcolor(goodTagIdx);
pdText = [{'None'};lc(:)];
colorMap = movieWindowHandles.colorMap;
textPos = def_textPos;
pd1Pos = def_pd1Pos;
pd2Pos = def_pd2Pos;
pdHandles = zeros(nGoodSpots,2);
% loop through spots and write name and two pulldowns
for iSpot = 1:nGoodSpots
    currentSpot = goodSpotsIdx(iSpot);

    % set positions
    % lines have a border of 0.15 at top and bottom and total height 1.7
    currentY = def_minY + (iSpot-1) * 2 + 0.15;
    [textPos(2),pd1Pos(2),pd2Pos(2)] = deal(currentY);


    % Look for fusions
    currentTagList = ...
        find(linklist(:,2) == linklist(currentSpot,2));
    if length(currentTagList) > 1
        % fusion spot. Use two colors for text, write 'fus'
        secondaryFusion = (linklist(currentTagList,3) == 4);

        % update goodSpotIdx, in case the good spot is the secondary
        % fusion
        changeIdx = goodSpotsIdx == currentTagList(secondaryFusion);
        if any(changeIdx)
            goodSpotsIdx(changeIdx) = currentTagList(~secondaryFusion);
        end

        th = uicontrol(reAssignH, 'Style', 'text', ...
            'Units','Characters', 'Position', textPos,...
            'FontName','Helvetica','HorizontalAlignment','center',...
            'FontWeight','bold',...
            'String','fusion',...
            'BackgroundColor',colorMap(currentTagList(~secondaryFusion),:),...
            'ForegroundColor',colorMap(currentTagList(secondaryFusion),:));
        pdHandles(iSpot,1) = uicontrol(reAssignH, 'Style','popupmenu',...
            'Units','Characters', ...
            'Position', pd1Pos, 'FontName','Helvetica','String',pdText,...
            'Value',find(currentTagList(~secondaryFusion)==goodTagIdx) + 1,...
            'Callback','LG_rag_tagPD_Callback(gcbo,[],guidata(gcbo));');
        pdHandles(iSpot,2) = uicontrol(reAssignH, 'Style','popupmenu',...
            'Units','Characters', ...
            'Position', pd2Pos, 'FontName','Helvetica','String',pdText,...
            'Value',find(currentTagList(secondaryFusion)==goodTagIdx) + 1,...
            'Callback','LG_rag_tagPD_Callback(gcbo,[],guidata(gcbo));');
    else
        % normal spot
        th = uicontrol(reAssignH, 'Style', 'text', ...
            'Units','Characters', 'Position', textPos,...
            'FontName','Helvetica','HorizontalAlignment','center',...
            'FontWeight','bold',...
            'String',idlistData.labelcolor(currentSpot),...
            'BackgroundColor',colorMap(currentSpot,:));
        pdHandles(iSpot,1) = uicontrol(reAssignH, 'Style','popupmenu',...
            'Units','Characters', ...
            'Position', pd1Pos, 'FontName','Helvetica','String',pdText,...
            'Value',find(currentTagList==goodTagIdx) + 1,...
            'Callback','LG_rag_tagPD_Callback(gcbo,[],guidata(gcbo));');
        pdHandles(iSpot,2) = uicontrol(reAssignH, 'Style','popupmenu',...
            'Units','Characters', ...
            'Position', pd2Pos, 'FontName','Helvetica','String',pdText,...
            'Value',1,...
            'Callback','LG_rag_tagPD_Callback(gcbo,[],guidata(gcbo));');
    end

end % loop spots

% add first line
currentY = def_minY + (nGoodSpots) * 2 + 0.15;
[textPos(2),pd1Pos(2),pd2Pos(2)] = deal(currentY);
th = uicontrol(reAssignH, 'Style', 'text', ...
    'Units','Characters', 'Position', textPos,...
    'FontName','Helvetica','HorizontalAlignment','center',...
    'FontWeight','normal',...
    'String',sprintf('t = %i',currentTime));
th = uicontrol(reAssignH, 'Style', 'text', ...
    'Units','Characters', 'Position', pd1Pos,...
    'FontName','Helvetica','HorizontalAlignment','center',...
    'FontWeight','normal',...
    'String','Tag 1');
th = uicontrol(reAssignH, 'Style', 'text', ...
    'Units','Characters', 'Position', pd2Pos,...
    'FontName','Helvetica','HorizontalAlignment','center',...
    'FontWeight','normal',...
    'String','Tag 2');


% --- updata handle structures

% assign default autput for reAssignGUI, remember handles of all pulldowns
handles.pdHandles = pdHandles;
pdValues = get(pdHandles,'Value');
pdValues = reshape([pdValues{:}],size(pdHandles));
handles.oldPdValues = pdValues;
handles.originalPdValues = pdValues;
handles.goodSpotsIdx = goodSpotsIdx;
handles.goodTagIdx = goodTagIdx;
else
    % assign pdValues, originalPdValues in case someone hits 'ok'
    handles.pdValues = 0;
    handles.originalPdValues = 0;
    
end % if ~noTags

% store currentTime
handles.currentTime = LG_getCurrentTime;

% Update handles structure
guidata(reAssignH, handles);
% Remember handle in movieWindowHandles
movieWindowHandles.otherWindows.LG_reAssignGUI = hObject;
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
% Remember position in naviHandles
naviHandles.positions.LG_reAssignGUI = newPosition;
guidata(naviHandles.LG_navigator, naviHandles);




% --- Outputs from this function are returned to the command line.
function varargout = LG_reAssignGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
varargout{1} = handles.output;
else
    varargout{1} = [];
end

% --- Executes on button press in LG_rag_cancel_pb.
function LG_rag_cancel_pb_Callback(hObject, eventdata, handles)

% don't ask. Just quit
close(handles.LG_reAssignGUI);


% --- Executes on button press in LG_rag_OK_pb.
function LG_rag_OK_pb_Callback(hObject, eventdata, handles)

% update idlist, no recalc
LG_rag_OK_Callback(handles,0);

function LG_rag_update_pb_Callback(hObject, eventdata, handles)

% close the gui and start over
close(handles.LG_reAssignGUI);
LG_reAssignGUI;


function LG_rag_OkAndRecalc_pb_Callback(hObject, eventdata, handles)

% update idlist, recalc
LG_rag_OK_Callback(handles,1);




% --- Executes on button press in LG_rag_futureFrames_RB.
function LG_rag_futureFrames_rb_Callback(hObject, eventdata, handles)
% hObject    handle to LG_rag_futureFrames_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn off other one if this one is turned on
myValue = get(hObject,'Value');
set(handles.LG_rag_thisFrame_rb,'Value',1-myValue)


% --- Executes on button press in LG_rag_thisFrame_RB.
function LG_rag_thisFrame_rb_Callback(hObject, eventdata, handles)
% hObject    handle to LG_rag_thisFrame_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn off other one if this one is turned on
myValue = get(hObject,'Value');
set(handles.LG_rag_futureFrames_rb,'Value',1-myValue)




