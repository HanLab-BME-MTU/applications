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
currentPath=get(handles.pathEdit,'String');
if isempty(currentPath)
    % No path specified yet
    set(handles.pathEdit,'String',fsmParam.main.path);
else
    if strcmp(fsmParam.main.path,get(handles.pathEdit,'String'))==0
    if ~isempty(fsmParam.main.path)
            infoPathString=['The loaded fsmParam.mat file contains stored path information [',fsmParam.main.path,'] which does not match the current selection [',currentPath,'].'];
            choice=questdlg(infoPathString, ...
                'User input requested', ...
                'Use current path','Use stored path','Use current path');
            switch choice
                case 'Use current path', set(handles.pathEdit,'String',currentPath);
                case 'Use stored path',  set(handles.pathEdit,'String',fsmParam.main.path);
                otherwise
                    error('Wrong selection');
            end
        else
            % No path specified yet
            set(handles.pathEdit,'String',fsmParam.main.path);
        end
    end
end

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
    
% Read parameter structure from handles.expPopup
fsmExpParam=get(handles.expPopup,'UserData');

% If the user selected 'Scale space', let's skip the whole settings check
if strcmp(fsmParam.main.label,'Scale space')
    fsmExpParam=[];
end

% If only the default fsmParam.mat has been loaded, there is no need for further controls
if ~isempty(fsmExpParam)
    
    % Check that the saved experiment number does not exceed the number of experiments in the database
    expNotValid=0;
    if fsmParam.main.noiseParam(7)>length(fsmExpParam)+1
        
        % Mark this experiment as not valid
        expNotValid=1;
        
    else
        
        % If the experiment number is valid, check that it is still pointing to the correct experiment
        if ~isfield(fsmParam.main,'label') % Back-compatibility
            fsmParam.main.label=[];
        end
        
        % This prevented a problem when the user attempts to start a new fsmGuiMain when it is already open
        %   (it is no longer strictly necessary)
        if fsmParam.main.noiseParam(7)==1
            expLabel='Select experiment';
        else
            expLabel=fsmExpParam(fsmParam.main.noiseParam(7)-1).label;
        end
        
        if expNotValid==0 & ~strcmp(fsmParam.main.label,expLabel)
            
            % Mark this experiment as not valid
            expNotValid=1;
            
        end
        
    end
    
    % If the experiment is valid, add it to the scroll-down menu, otherwise inform the user

    if expNotValid==0
        
        % Set correct experiment in scroll-down menu
        set(handles.expPopup,'Value',fsmParam.main.noiseParam(7));


        % Set correct experiment description
        if fsmParam.main.noiseParam(7)==1
            set(handles.textDescr,'String','Experiment description');
        else
            set(handles.textDescr,'String',fsmExpParam(fsmParam.main.noiseParam(7)-1).description);
        end

        % Check that that the noise parameters did not change
        if any((fsmParam.main.noiseParam(2:4)==fsmExpParam(fsmParam.main.noiseParam(7)-1).noiseParams)==0) | fsmParam.prep.gaussRatio~=fsmExpParam(fsmParam.main.noiseParam(7)-1).gaussRatio
            
            msg=['The parameters in the experiment database do not match those saved in the project. Which version do you want to use?'];
            choice=myQuestdlg(msg,'User input requested','Database (project parameters will be lost)','Project (database won''t be changed)','Database (project parameters will be lost)');
            if strcmp(choice,'Database (project parameters will be lost)')
                
                fsmParam.main.noiseParam(2:4)=fsmExpParam(fsmParam.main.noiseParam(7)-1).noiseParams;
                fsmParam.prep.gaussRatio=fsmExpParam(fsmParam.main.noiseParam(7)-1).gaussRatio;
                
            end

        end
        
        % If the noise parameters have been optimized, disable the quantile selection...
        if fsmExpParam(fsmParam.main.noiseParam(7)-1).quantile~=0
            set(handles.editZValue,'String',num2str(fsmExpParam(fsmParam.main.noiseParam(7)-1).quantile));
            fsmGuiUpdateConfidences(0);
        else
            
            % Otherwise make sure they are enabled
            set(handles.editZValue,'String','1.96');
            fsmGuiUpdateConfidences(1);
             
        end    
        
    else
    
        uiwait(msgbox('The experiment is no longer in the database (or old version of fsmParam.mat).','Error','modal'));
        set(handles.expPopup,'Value',1);
        set(handles.textDescr,'String','Experiment description');
        
    end

else

    if fsmParam.main.noiseParam(7)~=1
        
        error('Corrupted fsmParam.mat');
        
    else
    
        set(handles.expPopup,'Value',1);
        set(handles.textDescr,'String','Experiment description');

    end
    
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
    set(handles.textGauss,'Enable','on');
    set(handles.editGauss,'Enable','on');    
    set(handles.textDescr,'Enable','on');
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
    set(handles.textGauss,'Enable','off');
    set(handles.editGauss,'Enable','off');    
    set(handles.textDescr,'Enable','off');
end

% Gauss Ratio
set(handles.editGauss,'String',num2str(gaussRatio));

% Sigma (the following check assures back-compatibility)
if isfield(fsmParam.prep,'sigma')
    set(handles.editSigma,'String',num2str(fsmParam.prep.sigma));
else
    set(handles.editSigma,'String',num2str(1));    
end

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
    if fsmParam.track.tracker==1
        set(handles.editInfluence,'Enable','on');
    else
        set(handles.editInfluence,'Enable','off');
    end
    
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
    set(handles.editInfluence,'Enable','off');
end
% Default
set(handles.editThreshold,'String',num2str(fsmParam.track.threshold));
if isfield(fsmParam.track,'influence')
    set(handles.editInfluence,'String',num2str(fsmParam.track.influence));
else
    set(handles.editInfluence,'String',num2str(fsmParam.track.threshold));
end    
set(handles.checkGrid,'Value',fsmParam.track.grid);
if isfield(fsmParam.track,'corrLength')
    set(handles.editCorrLength,'String',num2str(fsmParam.track.corrLength));
else
    set(handles.editCorrLength,'String','33');
end
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

% Check for 'Scale space'
if fsmParam.prep.pstSpeckles==3 & strcmp(fsmParam.main.label,'Scale space')
    fsmParam.build.enable=0; % Make sure that only preprocessing and tracking are allowed if the user picked 'Scale space'
    set(handles.checkBuildModule,'Enable','off');
end

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

% Check for 'Scale space'
if fsmParam.prep.pstSpeckles==3 & strcmp(fsmParam.main.label,'Scale space')
    fsmParam.kin.enable=0; % Make sure that only preprocessing and tracking are allowed if the user picked 'Scale space'
    set(handles.textBleach,'Enable','off');
    set(handles.checkKinModule,'Enable','off');
    set(handles.bleachRadioOff,'Enable','off');
    set(handles.bleachRadio1x,'Enable','off');
    set(handles.bleachRadio2x,'Enable','off');
    set(handles.bleachRadio3x,'Enable','off');
end

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

% Check for 'Scale space'
if fsmParam.prep.pstSpeckles==3 & strcmp(fsmParam.main.label,'Scale space')
    fsmParam.disp.enable=0; % Make sure that only preprocessing and tracking are allowed if the user picked 'Scale space'
    set(handles.checkDispModule,'Enable','off');
end

% Enable/disable module
if fsmParam.disp.enable==1
    set(handles.checkDispModule,'Value',1);
else
    set(handles.checkDispModule,'Value',0);
end
