function fsmGuiWriteParameters(fsmParam,handles)
% fsmGuiWriteParameters sets the parameters read from fsmParam in the fsmGuiMain GUI 
%
% SYNOPSIS      fsmGuiWriteParameters(fsmParam,handles)
%
% INPUT         fsmParam :   general parameter structure for fsm
%                            Type help fsmGetParamDflts for more info on fsmParam
%               handles  :   fsmGuiMain GUI handles structure
%
% OUTPUT        NONE   
%                            
%
% DEPENDENCES              
%
% Aaron Ponti, November 1st, 2002

if nargin~=2
    error('Two input parameters expected');
end
    
% Constants
%confidenceProb=[1.15 1.29 1.45 1.645 1.96 2.58];
gaussRatio=fsmParam.prep.gaussRatio;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MAIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path
set(handles.pathEdit,'String',fsmParam.main.path);
% Image number
if fsmParam.specific.imageNumber~=0
    set(handles.numberEdit,'String',fsmParam.specific.imageNumber);
else
    set(handles.numberEdit,'String',fsmParam.main.imgN);
end
% Camera bit depth
set(handles.bitDepthEdit,'String',num2str(log2(fsmParam.main.normMax+1)));
% Noise parameters
set(handles.confOne,'Value',0);
set(handles.confTwo,'Value',0);
set(handles.confThree,'Value',0);
set(handles.confFour,'Value',0);
set(handles.confFive,'Value',0);
set(handles.confSix,'Value',0);
% Read z value from z-value edit box
zValue=get(handles.editZValue,'UserData');
switch fsmParam.main.noiseParam(6)
    case 1, set(handles.confOne,'Value',1); set(handles.editZValue,'String',num2str(zValue(1)));
    case 2, set(handles.confTwo,'Value',1); set(handles.editZValue,'String',num2str(zValue(2)));
    case 3, set(handles.confThree,'Value',1); set(handles.editZValue,'String',num2str(zValue(3)));
    case 4, set(handles.confFour,'Value',1); set(handles.editZValue,'String',num2str(zValue(4)));
    case 5, set(handles.confFive,'Value',1); set(handles.editZValue,'String',num2str(zValue(5)));
    case 6, set(handles.confSix,'Value',1); set(handles.editZValue,'String',num2str(zValue(6)));
    otherwise, set(handles.editZValue,'String',fsmParam.main.noiseParam(5)); % User-defined quantile
end
% Set correct experiment in scroll-down menu
set(handles.expPopup,'Value',fsmParam.main.noiseParam(7));

% Read parameter structure from handles.expPopup
fsmExpParam=get(handles.expPopup,'UserData');

% Set correct experiment description
if fsmParam.main.noiseParam(7)==1
    set(handles.textDescr,'String','Experiment description');
else
    set(handles.textDescr,'String',fsmExpParam(fsmParam.main.noiseParam(7)-1).description);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PREPROCESSING MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable/disable module
if fsmParam.prep.enable==1
    set(handles.checkPrepModule,'Value',1);
else
    set(handles.checkPrepModule,'Value',0);
end
% Default values
switch fsmParam.prep.pstSpeckles
case 1
    set(handles.primaryRadio,'Value',1);
    set(handles.tertiaryRadio,'Value',0);
    set(handles.scaleRadio,'Value',0);
    set(handles.percText,'Enable','off');
    set(handles.percEdit,'Enable','off');
    set(handles.orderText,'Enable','off');
    set(handles.orderEdit,'Enable','off');
    set(handles.sigText,'Enable','off');
    set(handles.sigEdit,'Enable','off');
case 2
    set(handles.primaryRadio,'Value',0);
    set(handles.tertiaryRadio,'Value',1);
    set(handles.scaleRadio,'Value',0);
    set(handles.percText,'Enable','on');
    set(handles.percEdit,'Enable','on');
    set(handles.orderText,'Enable','on');
    set(handles.orderEdit,'Enable','on');
    set(handles.sigText,'Enable','off');
    set(handles.sigEdit,'Enable','off');
case 3
    set(handles.primaryRadio,'Value',0);
    set(handles.tertiaryRadio,'Value',0);
    set(handles.scaleRadio,'Value',1);
    set(handles.percText,'Enable','off');
    set(handles.percEdit,'Enable','off');
    set(handles.orderText,'Enable','off');
    set(handles.orderEdit,'Enable','off');
    set(handles.sigText,'Enable','on');
    set(handles.sigEdit,'Enable','on');
otherwise
    error('Value for PST out of range.');
end

set(handles.TriangCheck,'Value',fsmParam.prep.enhTriang);
set(handles.autoPolCheck,'Value',fsmParam.prep.autoPolygon);
set(handles.drawROICheck,'Value',fsmParam.prep.drawROI);

set(handles.orderEdit,'String',num2str(fsmParam.prep.paramSpeckles(1))); % Sets the order for 'higher-order speckles'
set(handles.percEdit,'String',num2str(fsmParam.prep.paramSpeckles(2)));   % Sets the percentage for 'higher-order speckles'
set(handles.sigEdit,'String',num2str(fsmParam.prep.paramSpeckles(3)));   % Sets the sigma for 'scale space speckles'

if fsmParam.prep.enable==1
    set(handles.TriangCheck,'Enable','on');
    set(handles.autoPolCheck,'Enable','on');
    set(handles.textDel,'Enable','on');
    set(handles.textCameraCalPar,'Enable','on');    
    set(handles.expPopup,'Enable','on');
    set(handles.primaryRadio,'Enable','on');
    set(handles.tertiaryRadio,'Enable','on');
    set(handles.textAdvancedPrep,'Enable','on');
    set(handles.expPopup,'Enable','on');
    set(handles.textGauss,'Enable','on');
    set(handles.editGauss,'Enable','on');
    set(handles.textExplanation,'Enable','on');
    set(handles.textExplanation2,'Enable','on');
else
    set(handles.TriangCheck,'Enable','off');
    set(handles.autoPolCheck,'Enable','off');
    set(handles.textDel,'Enable','off');
    set(handles.textCameraCalPar,'Enable','off');
    set(handles.expPopup,'Enable','on');
    set(handles.primaryRadio,'Enable','off');
    set(handles.tertiaryRadio,'Enable','off')
    set(handles.textAdvancedPrep,'Enable','off');
    set(handles.expPopup,'Enable','off');
    set(handles.textGauss,'Enable','off');
    set(handles.editGauss,'Enable','off');
    set(handles.textExplanation,'Enable','off');
    set(handles.textExplanation2,'Enable','off');  
end

% Gauss Ratio
set(handles.editGauss,'String',num2str(gaussRatio));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% TRACKING MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable/disable module
if fsmParam.track.enable==1
    set(handles.checkTrackModule,'Value',1);
    set(handles.checkEnhTrack,'Enable','on');
    set(handles.radioTrackBrownian,'Enable','on');
    set(handles.radioTrackBrownian,'Enable','on');
    set(handles.radioTrackBrownian,'Enable','on');
    set(handles.radioEnhTrackBrownian,'Enable','on');
    if fsmParam.track.enhanced==1
        set(handles.checkGrid,'Enable','on');
    else
        set(handles.checkGrid,'Enable','off');
    end               
else
    set(handles.checkTrackModule,'Value',0);
    set(handles.radioTrackBrownian,'Enable','off');
    set(handles.radioEnhTrackBrownian,'Enable','off');
    set(handles.radioTrackFlow,'Enable','off');
    set(handles.textThreshold,'Enable','off');
    set(handles.editThreshold,'Enable','off');
    set(handles.checkEnhTrack,'Enable','off');
    set(handles.checkGrid,'Enable','off');
end
% Default
set(handles.editThreshold,'String',fsmParam.track.threshold);
set(handles.checkEnhTrack,'Value',fsmParam.track.enhanced);
set(handles.checkGrid,'Value',fsmParam.track.grid);
switch fsmParam.track.tracker
    case 1
        set(handles.radioTrackBrownian,'Value',1);
        set(handles.radioEnhTrackBrownian,'Value',0);
        set(handles.radioTrackFlow,'Value',0);
    case 2
        set(handles.radioTrackBrownian,'Value',0);
        set(handles.radioEnhTrackBrownian,'Value',1);
        set(handles.radioTrackFlow,'Value',0);
    case 3
        set(handles.radioTrackBrownian,'Value',0);
        set(handles.radioEnhTrackBrownian,'Value',0);
        set(handles.radioTrackFlow,'Value',1);
otherwise
        error('Please specify an existing tracker');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BUILDER MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable/disable module
if fsmParam.build.enable==1
    set(handles.checkBuildModule,'Value',1);
else
    set(handles.checkBuildModule,'Value',0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% KINETIC ANALYSIS MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable/disable module
if fsmParam.kin.enable==1
    set(handles.checkKinModule,'Value',1);
else
    set(handles.checkKinModule,'Value',0);
end
% Set bleaching reduction to off
set(handles.bleachRadioOff,'Value',0);
set(handles.bleachRadio1x,'Value',0);
set(handles.bleachRadio2x,'Value',0);
set(handles.bleachRadio3x,'Value',0);
switch fsmParam.kin.bleachRed
    case 0,        set(handles.bleachRadioOff,'Value',1);
    case 7.25e-5,  set(handles.bleachRadio1x,'Value',1);
    case 1.45e-4,  set(handles.bleachRadio2x,'Value',1);
    case 2.175e-4, set(handles.bleachRadio3x,'Value',1);
    otherwise error('Wrong bleaching value');
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RESULT DISPLAY MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable/disable module
if fsmParam.disp.enable==1
    set(handles.checkDispModule,'Value',1);
else
    set(handles.checkDispModule,'Value',0);
end
