function fsmParam=fsmGetParamDflts
% fsmGetParamDflts returns the default fsmParam structure to be used with fsmMain
%
% SYNOPSIS       fsmParam=fsmGetParamDflts
%
% INPUT          None
%
% OUTPUT         fsmParam: default fsmParam structure for fsmMain
%
% DEPENDENCIES   fsmGetParamDflts uses { nothing }
%                fsmGetParamDflts is used by { fsmMain }
%
% Aaron Ponti, October 2nd, 2002

if nargin~=0
    error('fsmgetParamDflts accepts no input parameter');
end

% Create structure and substructures
fsmParam=struct('main',0,'prep',0,'track',0,'build',0,'kin',0,'disp',0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MAIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path
fsmParam.main.path='';                 % Work and result path (select from GUI)
% Image path
fsmParam.main.imagePath='';            % Image path (defined at the project level)
% Image number
fsmParam.main.imgN=0;                  % Number of image to be processed from the stack [0 = all]
% Camera bit depth
fsmParam.main.normMin=0;               % Lower intensity bound for intensity normalization (gray values) [0]
fsmParam.main.normMax=2^14-1;          % Upper intensity bound for intensity normalization (gray values) [2^14-1]
% Noise paramaters
fsmParam.main.noiseParam=[0 0 0 0 0 5 1];% Noise parameters from the camera calibration (select from GUI)
                                         % The sixth position defines which confidence was chosen
                                         % The seventh position which noise parameter wae set
% To prevent pointing to wrong experiments in the database, store the label of the experiment, too
fsmParam.main.label='';
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PREPROCESSING MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable module
fsmParam.prep.enable=1;                    % Activates the PREPROCESSING module [ 0 : off, 1 : on ]
% Secondary speckles
fsmParam.prep.pstSpeckles=1;               % Sets the type of speckles to be detected
fsmParam.prep.paramSpeckles=[3 0.01 1.06]; % Sets [order percentage] for 'higher-order speckles' and sigma for 'scale space spaeckles'

% Enhanced triangulation
fsmParam.prep.enhTriang=0;          % Delaunay triangulation can be made more stable by enhancing it 
%                                    [ 0 : no | 1 : yes (slower)]
% Auto polygon
fsmParam.prep.autoPolygon=0;        % Automatic analisys of the image to extract cell boundaries  [ 0 : off, 1 : on ]
% The user draws a ROI
fsmParam.prep.drawROI=0;            % The user can draw a ROI on the image to restrict the analysis

fsmParam.prep.gaussRatio=3.54;      % Gauss ratio
fsmParam.prep.sigma=1;              % Sigma for low-pass filtering
fsmParam.prep.uptodate=0;           % Up-to-date flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% TRACKING MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable module
fsmParam.track.enable=1;            % Activates the TRACKER module [ 0 : off | 1 : on ]
% Threshold
fsmParam.track.threshold=3;         % Defines the search radius for the tracker (pixels)
fsmParam.track.influence=3;         % Defines the radius of influence for the BMTNN tracker (pixels)
% Tracker
fsmParam.track.tracker=1;           % Specifies the tracker to be selected [1: BMTNN | 2: BMTG | 3: 3FT ]
% Hierarchical tracking
fsmParam.track.enhanced=0;          % When this is flag is set to one, the tracking is performed twice
                                    % The interpolated vector field obtained in the first tracking is
                                    % used to refine the tracking for the second round
% Hierarchical tracking on a grid
fsmParam.track.grid=0;              % Specifies whether the interpolation of the vector field has to be
                                    % be performed onto a regular grid or at the original
                                    % vector positions
% Correlation length for tracking
fsmParam.track.corrLength=33;       % Correlation length in pixels

fsmParam.track.uptodate=0;          % Up-to-date flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BUILDER MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable module
fsmParam.build.enable=1;            % Activates the BUILDER module [ 0 : off | 1 : on ]

fsmParam.build.uptodate=0;          % Up-to-date flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% KINETIC ANALYSIS MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable module
fsmParam.kin.enable=1;              % Activates the KINETIC ANALYSIS module [ 0 : off | 1 : on ]
% Bleaching reduction
fsmParam.kin.bleachRed=0;           % Bleaching reduction filters out weak scores (in the bleaching range) [ 0 | 1 | 2 | 3]
                                    % [ 0        : no selection
                                    %   7.25e-5 : ~0.5*mean(bleaching_scores in 161fix)
                                    %   1.45e-4  : ~1*mean(bleaching_scores in 161fix)
                                    %  2.175e-4 : ~1.5*mean(bleaching_scores in 161fix)

fsmParam.kin.uptodate=0;          % Up-to-date flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RESULT DISPLAY MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable module
fsmParam.disp.enable=1;             % Activates the  RESULT DISPLAY module [ 0 : off | 1 : on ]

fsmParam.disp.uptodate=0;           % Up-to-date flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% EXPERIMENT SPECIFIC - this entries are generated by the software during
%    runtime
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsmParam.specific.imgSize=[0 0];    % Image size [y x]
fsmParam.specific.imageNumber=0;    % Number of images to be treated (the actual number)
fsmParam.specific.formString='';    % String to correctly format the numeric suffix of saved files
fsmParam.specific.fileList='';      % All file names written in a (char) matrix; eah row contains 
                                    % a file name with full path
fsmParam.specific.intCorrFactors=0; % Intensity correction factors
fsmParam.specific.firstIndex=0;     % Index of the first image
fsmParam.specific.lastIndex=0;      % Index of the last image
fsmParam.specific.lastRun='';       % Date and time of the last run