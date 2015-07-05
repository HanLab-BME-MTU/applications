function varargout = Gui_AnalyzeMovie(varargin)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

% GUI_ANALYZEMOVIE M-file for Gui_AnalyzeMovie.fig
%      GUI_ANALYZEMOVIE, by itself, creates a new GUI_ANALYZEMOVIE or raises the existing
%      singleton*.
%
%      H = GUI_ANALYZEMOVIE returns the handle to a new GUI_ANALYZEMOVIE or the handle to
%      the existing singleton*.
%
%      GUI_ANALYZEMOVIE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ANALYZEMOVIE.M with the given input arguments.
%
%      GUI_ANALYZEMOVIE('Property','Value',...) creates a new GUI_ANALYZEMOVIE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gui_AnalyzeMovie_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_AnalyzeMovie_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Gui_AnalyzeMovie

% Last Modified by GUIDE v2.5 25-Mar-2002 12:19:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_AnalyzeMovie_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_AnalyzeMovie_OutputFcn, ...
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


% --- Executes just before Gui_AnalyzeMovie is made visible.
function Gui_AnalyzeMovie_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gui_AnalyzeMovie (see VARARGIN)
% Choose default command line output for Gui_AnalyzeMovie
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    % UIWAIT makes Gui_AnalyzeMovie wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    imagesFolder = 'C:\Documents and Settings\godina\Desktop\Skeletons\';
    analyzeSkeletonMovie
    
% --- Outputs from this function are returned to the command line.
function varargout = Gui_AnalyzeMovie_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
