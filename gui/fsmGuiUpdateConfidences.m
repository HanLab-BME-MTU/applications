function fsmGuiUpdateConfidences(enable)
% fsmGuiUpdateConfidences reads the quantile from fsmGuiMain and updates the confidence radio buttons
%
% SYNOPSIS      fsmGuiUpdateCondifences
%
% INPUT         handles : handles of the fsmGuiMain GUI (set equal to [] to have the function search for them)
%               enable  : [ 0 | 1 ] 0 : disables all confidence radio buttons and the quantile edit field
%                                   1 : enables all confidence radio buttons and the quantile edit field 
%                         This is used to easily turn on/off the selection of confidence probabilities for 
%                         non-optimized/optimized experimental settings.
%
% OUTPUT        none
%
% DEPENDENCES   fsmGuiUpdateConfidences uses { }
%               fsmGuiUpdateCondifences is used by { fsmGuiMain
%                                                    fsmGuiWriteParameters }
%
% Aaron Ponti, May 14th, 2004

% Find fsmGuiMain
h=findall(0,'Tag','fsmGuiMain');

if isempty(h)
    disp('fsmGuiMain not running');
    return
end

% Get all handles
handles=guidata(h);

% Set all radio buttons to 0
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);

% Read the list of quantile presets
zValues=get(handles.editZValue,'UserData');

% Read the entered quantile
z=str2num(get(handles.editZValue,'String'));

% Check whether the entered quantile is one of the presets
indx=find(zValues==z);
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);
if ~isempty(indx)
    switch(indx)
        case 1, set(handles.confOne,'Value',1);
        case 2, set(handles.confTwo,'Value',1);
        case 3, set(handles.confThree,'Value',1);
        case 4, set(handles.confFour,'Value',1);
        case 5, set(handles.confFive,'Value',1);
        case 6, set(handles.confSix,'Value',1);
        otherwise, error('Somehow the entered value does match one of the presets but still we can find which');
    end
end

% Turn on/off all radio buttons and edit fiels
if enable==0
    set(handles.confOne,'Enable','off');
    set(handles.confTwo,'Enable','off');
    set(handles.confThree,'Enable','off');
    set(handles.confFour,'Enable','off');
    set(handles.confFive,'Enable','off');
    set(handles.confSix,'Enable','off');
    set(handles.editZValue,'Enable','off');
    set(handles.textConfidence,'Enable','off');
    set(handles.textZValue,'Enable','off');
else
    set(handles.confOne,'Enable','on');
    set(handles.confTwo,'Enable','on');
    set(handles.confThree,'Enable','on');
    set(handles.confFour,'Enable','on');
    set(handles.confFive,'Enable','on');
    set(handles.confSix,'Enable','on');
    set(handles.editZValue,'Enable','on');
    set(handles.textConfidence,'Enable','on');
    set(handles.textZValue,'Enable','on');
end
