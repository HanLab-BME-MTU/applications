function [fsmParam,status]=fsmGuiReadParameters(handles,fsmParam)
% fsmGuiReadParameters reads parameters from fsm user interface 
%
% SYNOPSIS      [fsmParam,status]=fsmGuiReadParameters(fsmParam,handles)
%
% INPUT         handles  :   structure containing the handles of the figure fsmGuiMain
%               fsmParam :   (optional) general parameter structure for fsm
%                            If fsmParam is not passed as an input, the structure is 
%                            initialized and filled with the values read from the min GUI.
%                            Otherwise, the structure is completed. 
%                            Type help fsmGetParamDflts for more info on 
%                            fsmParam
%
% OUTPUT        fsmParam :   generated/augmented fsmParam
%               status   :   1 if fsmParam was successfully filled/generated,
%                            and 0 otherwise.
%
% DEPENDENCES   
%               
%
% Aaron Ponti, October 2nd, 2002

if nargin==0
    % Initialize structure and substructures
    fsmParam=struct('main',0,'prep',0,'track',0,'build',0,'kin',0,'disp',0);
end

% Set status
status=0;

% Constants
%confidenceProb=[1.15 1.29 1.45 1.645 1.96 2.58];
gaussRatio=str2num(get(handles.editGauss,'string'));

% Modules
fsmParam.prep.enable=get(handles.checkPrepModule,'Value');
fsmParam.track.enable=get(handles.checkTrackModule,'Value');
fsmParam.build.enable=get(handles.checkBuildModule,'Value');
fsmParam.kin.enable=get(handles.checkKinModule,'Value');
fsmParam.disp.enable=get(handles.checkDispModule,'Value');

% Path
path=get(handles.pathEdit,'String');
if isempty(path)
    uiwait(msgbox('Please setup a project in fsmCenter.','Warning','warn'));
    return;
else
    fsmParam.main.path=path;
end

% Image path
fsmParam.main.imagePath=get(handles.textImage,'String');

% Number of images
fsmParam.main.imgN=str2num(get(handles.numberEdit,'String'));

% Auto cell boundaries
fsmParam.prep.autoPolygon=get(handles.autoPolCheck,'Value');

% Special Bit depth setting for Auto cell boundaries
val                 = get(handles.edgeBitDepth,'Value');
string_list         = get(handles.edgeBitDepth,'String');
bit_depth           = string_list{val};
fsmParam.prep.edgeBitDepth = bit_depth;

% User draws a ROI
drawROI=get(handles.drawROICheck,'Value');
loadROI=get(handles.loadROICheck,'Value');
if drawROI==0 & loadROI==0
    fsmParam.prep.drawROI=0;    % No drawing, no loading
elseif drawROI==1 & loadROI==0
    fsmParam.prep.drawROI=1;    % Drawing, no loading
elseif drawROI==0 & loadROI==1
    fsmParam.prep.drawROI=2;    % No drawing, loading
else
    error('This combination of selections for load/draw ROI should be impossible.');
end

% Gauss ratio
fsmParam.prep.gaussRatio=str2num(get(handles.editGauss,'String'));

% Sigma
fsmParam.prep.sigma=str2num(get(handles.editSigma,'String'));

% Noise parameters
indx=find([get(handles.confOne,'Value') get(handles.confTwo,'Value') get(handles.confThree,'Value') get(handles.confFour,'Value') get(handles.confFive,'Value') get(handles.confSix,'Value')]);
fsmParam.main.noiseParam(1,1)=str2num(get(handles.editZValue,'String'))/gaussRatio;
fsmParam.main.noiseParam(1,5)=str2num(get(handles.editZValue,'String'));
if ~isempty(indx)
    fsmParam.main.noiseParam(1,6)=indx; % Stores the index of the radio button to be set to 1
else
    fsmParam.main.noiseParam(1,6)=0;
end

% Read parameter structure from handles.expPopup
fsmExpParam=get(handles.expPopup,'UserData');

% Add value of the popup window
fsmParam.main.noiseParam(1,7)=get(handles.expPopup,'Value');

% Normalization boundaries
fsmParam.main.normMin=0;
fsmParam.main.normMax=2^(str2num(get(handles.bitDepthEdit,'String')))-1;

% Bleaching threshold for selecting scores
bleachValues=[0 7.25e-5 1.45e-4 2.175e-4]; % 1x corresponds to the mean of the negative bleaching scores
indx=find([get(handles.bleachRadioOff,'Value') get(handles.bleachRadio1x,'Value') get(handles.bleachRadio2x,'Value') get(handles.bleachRadio3x,'Value')]);
fsmParam.kin.bleachRed=bleachValues(indx);


% Subpixel accurate localization of speckles
fsmParam.prep.subpixel = get(handles.subpixel,'Value');

% Secondary speckles
fsmParam.prep.pstSpeckles=find([get(handles.primaryRadio,'Value') get(handles.tertiaryRadio,'Value') get(handles.scaleRadio,'Value')]);
switch fsmParam.prep.pstSpeckles
    case 1
        % Only primary speckles
        fsmParam.prep.paramSpeckles(1)=1; % Sets the order for 'higher-order speckles'
        fsmParam.prep.paramSpeckles(2)=str2num(get(handles.percEdit,'String')); % Percentage
        fsmParam.prep.paramSpeckles(3)=str2num(get(handles.sigEdit,'String')); % Sets the sigma for 'scale space speckles'
    case 2
        % Higher-order speckles
        fsmParam.prep.paramSpeckles(1)=str2num(get(handles.orderEdit,'String'));
        fsmParam.prep.paramSpeckles(2)=str2num(get(handles.percEdit,'String')); % Percentage
        fsmParam.prep.paramSpeckles(3)=str2num(get(handles.sigEdit,'String'));
    case 3
        % Scale space speckles
        fsmParam.prep.paramSpeckles(1)=str2num(get(handles.orderEdit,'String'));
        fsmParam.prep.paramSpeckles(2)=str2num(get(handles.percEdit,'String')); % Percentage
        fsmParam.prep.paramSpeckles(3)=str2num(get(handles.sigEdit,'String'));
    otherwise
        error('The expected value for fsmParam.prep.pstSpeckles is either 1 or 2 or 3');
end

% Assign noise parameters
exp=get(handles.expPopup,'Value');
if exp~=1
    fsmParam.main.noiseParam(1,2:4)=fsmExpParam(exp-1).noiseParams;
    % Store experiment label
    fsmParam.main.label=fsmExpParam(exp-1).label;
else
    if fsmParam.prep.pstSpeckles(1)~=3
        uiwait(msgbox('Please select an experiment.','Warning','warn'));
        return;
    else
        % Store experiment label
        fsmParam.main.label='Scale space'; % Set it to 'Scale space' so that there won't be any check for setting the next time fsmParam is loaded
        % Make sure that fsmParam.main.noiseParam(7) does not point to any experiment
        fsmParam.main.noiseParam(7)=1;
        set(handles.expPopup,'Value',1);
    end
end


% Enhanced triangulation
fsmParam.prep.enhTriang=get(handles.TriangCheck,'Value');

% Tracker
indx=get(handles.popupTracker,'Value');
fsmParam.track.tracker=indx;
fsmParam.track.enhanced=get(handles.checkEnhTrack,'Value');
fsmParam.track.grid=get(handles.checkGrid,'Value');

% Tracker's search radius
fsmParam.track.threshold=str2num(get(handles.editThreshold,'String'));
fsmParam.track.influence=str2num(get(handles.editInfluence,'String'));
fsmParam.track.initRadius=str2num(get(handles.editInitRadius,'String'));

% Tracker's correlation length
fsmParam.track.corrLength=str2num(get(handles.editCorrLength,'String'));

% Tracker's initializer
if get(handles.checkTrackInit,'Value')==0
    fsmParam.track.init=0;
    fsmParam.track.initPath='';
else
    fsmParam.track.init=get(handles.popupTrackInit,'Value');
    switch fsmParam.track.init
        case 1,
            % Vectors provided by imKymoAnalysis - saved into /corr/flow
            fsmParam.track.initPath=[fsmParam.project.path,filesep,fsmParam.project.corr,filesep,'flow',filesep];
          
        case 2,
            % Vectors provided by TFT - saved into /tack/flow
            fsmParam.track.initPath=[fsmParam.project.path,filesep,fsmParam.project.tack,filesep,'flow',filesep];
            
        otherwise
            
            error('Only tracker initialization by imKymoAnalysis or TFT is supported right now.');
    end
end

% Set status to 1
status=1;
