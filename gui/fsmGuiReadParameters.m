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
    uiwait(msgbox('Please choose a work directory.','Warning','warn'));
    return;
else
    fsmParam.main.path=path;
end

% Number of images
fsmParam.main.imgN=str2num(get(handles.numberEdit,'String'));

% Auto cell boundaries
fsmParam.prep.autoPolygon=get(handles.autoPolCheck,'Value');

% User draws a ROI
fsmParam.prep.drawROI=get(handles.drawROICheck,'Value');

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

% Assign noise parameters
exp=get(handles.expPopup,'Value');
if exp~=1
    fsmParam.main.noiseParam(1,2:4)=fsmExpParam(exp-1).noiseParams;
else
    uiwait(msgbox('Please select an experiment.','Warning','warn'));
    return;
end

% Add value of the popup window
fsmParam.main.noiseParam(1,7)=get(handles.expPopup,'Value');

% Store experiment label
fsmParam.main.label=fsmExpParam(exp-1).label;

% Normalization boundaries
fsmParam.main.normMin=0;
fsmParam.main.normMax=2^(str2num(get(handles.bitDepthEdit,'String')))-1;

% Bleaching threshold for selecting scores
bleachValues=[0 7.25e-5 1.45e-4 2.175e-4]; % 1x corresponds to the mean of the negative bleaching scores
indx=find([get(handles.bleachRadioOff,'Value') get(handles.bleachRadio1x,'Value') get(handles.bleachRadio2x,'Value') get(handles.bleachRadio3x,'Value')]);
fsmParam.kin.bleachRed=bleachValues(indx);

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

% Enhanced triangulation
fsmParam.prep.enhTriang=get(handles.TriangCheck,'Value');

% Tracker
indx=find([get(handles.radioTrackBrownian,'Value') get(handles.radioEnhTrackBrownian,'Value') get(handles.radioTrackFlow,'Value')]);
fsmParam.track.tracker=indx;
fsmParam.track.enhanced=get(handles.checkEnhTrack,'Value');
fsmParam.track.grid=get(handles.checkGrid,'Value');

% Tracker's search radius
fsmParam.track.threshold=str2num(get(handles.editThreshold,'String'));
fsmParam.track.influence=str2num(get(handles.editInfluence,'String'));

% Set status to 1
status=1;
