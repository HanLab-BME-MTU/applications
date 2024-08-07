%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imaging parameters: (acquired from 'lastProjSettings.mat')
% actPixelSize: pixel size of actin channel in nm.
% frameInterval: frame interval in seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% trackMethod: Specify the source of flow data. 
%    Possible values:
%    'corr' : correlation tracking
%    'speck': specktackle particle tracking in FSM.
%
% isFieldBndFixed: For dynamic force reconstruction, specify whether the field boundary 
%                  is fixed or not. Possible values: 'yes'(default) or 'no'. When it is 'yes', we
%                  only need to calculate the forward map for each basis function in the domain once
%                  for all time points.
% isDataPosFixed : For dynamic force reconstruction, specify whether the position of data points in
%                  the field is fixed. Possible values: 'yes' or 'no' (default). This variable only
%                  makes sense when 'isFieldBndFixed' is 'yes'. When it is 'yes', we only need to
%                  calculate the forward operator matrix 'A' once for all time points.
%
% YModulVariation: Specify whether homogeneous ('homo' default) or inhomogeneous ('inhomo' or 
%                  'inhomoUser' (user defined)) Young's modulus is used. When it is 'inhomo', 
%                  Young's modulus is proportional to image intensity.
%
% The following are parameters used for filtering and interpolation of raw data to grid points. 
%    corLen   : Used for a faithful interpolation with minimum filtering. Unit, pixels.
%    sCorLen  : Used to calculate smoothly interpolated field. Unit, pixels.
%    edgCorLen: For calculating the displacements on the boundary.
%
% gridDx: The grid size (in pixels) in the x-direction.
% gridDy: The grid size (in pixels) in the y-direction.
%
% dataSite: Specify the data sites where the forward propagated displacements and the
%           tracked displacements are compared. 
%    Possible Values:
%    'grid'      : On grid points.
%    'original'  : If 'trackMethod' is 'speck', it means the speckles in the 
%                  starting frame of each time point. If 'trackMethod' is
%                  'corr', it is whatever tracked points in 'flowTrack' data.
%    'everyOther': Sometime the origianlly tracked points are too dense.
%                  To reduce computation, you can specify to use every
%                  other points in the original.
%
% dataToUse: Specify what data to use. 
%    Possible Values : 
%    'interp' : Displacements interpolated from the original (raw) displacements
%               vectors with very small space correlation length: 'corLen'.
%    'smooth' : The filtered displacement vectors with relatively big correlation 
%               length: 'sCorLen'.
%    'simul'  : The displacement computed from the modified body force
%               by the script 'retroFlowSimulate'. This is for the debugging and
%               justification of the force identification program.
%
% bfDisplaySite: Specify where the identified body force is displayed. 
%    Possible Values:
%    'grid'      : On grid points.
%    'original'  : The original raw data points. For data from speckTackle, 
%                  these are the speckles in the starting frame of each time 
%                  step. For data from correlation tracking, These are the
%                  tracked points for the corresponding time step.
%    'everyOther': Sometime the origianlly tracked points are too dense.
%                  To increase the visibility, you can specify to use every
%                  other point in the original.
%    'dataSite'  : The same as the data points used in the reconstruction.
%
% The following two parameters specify if the forward operator has already been computed.
%    bfFwdOpComputed: For domain force.
%    tfFwdOpComputed: For boundary force.
%
%    Possible Values : 
%    'fem'  : The femlab solution for each basis function has been calculated.
%    'A'    : The matrix representation of the forward operator has been 
%             calculated.
%    'all'  : Everything is done.
%    'none' : Nothing has been calculated.
%
% forceToIdentify: Specify what force is to be identified.
%    Possible Values : 
%    'bf' : Body force.
%    'tf' : Boundary traction force.

% Specify the threshold angles between the force and the displacements that are
% used in the cone rule to determine if the force is mainly the contraction force or the
% adhesion force.
%    mcfAngle: For contraction force.
%    adfAngle: For boundary force.
%
% fixMcfProj: Specify whether the projected direction of contraction is fixed at
%             'mcfAngle'.
%
% edgBrkDistTF: The distance between the breaks of each edge that are used to build the
%               B-spline interpolation of the boundary traction force.
% bspOrderTF  : The order of spline for boundary force.
% edgDispIntvLen: The interval length between the points on each edge used
%                to show result.
% edgToDisplay: Edges to display boundary force.
%
% Specify the regularization parameter to be used in solving the least-square problem.
%    bfSigma: For domain force.
%    tfSimga: For boundary force.
%
% testTimeStep: To test the identification program, pick one time step and use the identified
%               force as input to 'retroFlowSimulate' to get simulated data.
%
% numBSolsPerFile: Because of memory issue, we save the solutions for basis functions in
%                  a number of files. Each file contains 'numBSolsPerFile'. It should be an
%                  even number.
% overwriteBSolBF: 'yes' (default) or 'no'. Whether overwrite previously calculated and saved basis
%                  solution for reconstructing body force.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters for displaying vector field, image and force color map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the scale to be used by 'quiver' to display the vector.
%    bfScale  : for the body force.
%    tfScale  : for the boundary force.
%    dispScale: for the displacement.
%
% imgChannel: Specify which image channel to be used for color map overlay.
% relDispImgFrmNo: Relative frame number of image for display.
%
% smDispThreshold: Specify a threshold that is used to determine if the displacement is too
%                  small to show the directionality.
%
% smForceThreshold: Small force threshold in terms of a percentage of the maximum force.
%
% Specify whether draw the boundary of contraction and adhesion region in 
% the force color map and what color to use.
%    markMcfBnd : Whether mark the boundary of contraction force. Default: 'yes';
%    markAdfBnd : Whether mark the boundary of adhesion resistance force. Default: 'yes';
%    mcfBndColor: Color for contraction force boundary. Default: 'w';
%    adfBndColor: Color for adhesion force boundary. Default: 'r';
%    showFlowVec: Whether flow vector is shown. Default: 'yes';
%    showBdfVec : Whether boundary force vector is shown. Default: 'yes';
%    showMcfVec : Whether contraction force vector is shown. Default: 'yes';
%    showAdfVec : Whehter adhesion force vector is shown. Default: 'yes';
%
% showMixZone: Specify whether show color map of separated forces in the mixed zone of 
%              'mcfMap' and 'adfMap'.
%
% markMixZone: Specify whether mark (with dotts and boundary) the mixed zone of in the 
%              separated contraction force and adhesion force map.
% markCellEdge: Specify whether cell edges are displayed. It only works if
%               edges have been detected and are saved in the corresponding
%               edge directory.
% cellEdgeColor: Specify the color of the displayed cell edge.
%
% mixfBndColor: Color of mixed region boundary.
%
%imgIRange: Image intensity range. To correct for the bright spot in some image. 
%           Default is [0 1] (full range).
%
% Force color map range.
%   bdfColorDispRange: Domain force. 
%   mcfColorDispRange: Contraction force.
%   adfColorDispRange: Adhesion force.
%   spdColorDispRange: Speed.
%
% showBdfCBar: Specify whether show color bar for the speed map.
% showSpdCBar: Specify whether show color bar for the speed map.
%
% showYModCBar: Specify whether show color bar for the Young's modulus color map.
%
% spdUnit: Specify speed unit.
%    Possible Value:
%    'pixelPerFrame', 'pixelPerMin', 'umPerMin' and 'nmPerMin'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter for retrograde flow simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ForceToSimu: Specify what forces are simulated.
%    Possible values:
%    'mcfOnly'  : Simulate only contraction force.
%    'mcfAndBnd': Contraction and boundary force.
%    'mcfAndAdf': Contraction and adhesion.
%    'adfOnly'  : Simulate only adhesion.
%    'adfAndBnd': Adhesion and boundary.
%    'bndOnly'  : Simulate only boundary force.
%    'all'      : All three types of forces: contraction, adhesion and
%                 boundary.
% simMCFMagFactor : The magnification factor over the recovered myosin 
%                   contraction force used in the simulation.
% simADCMagFactor : The magnification factor over the recovered adhesion
%                   coefficient used in the simulation.
% simAbsNoiseLevel: Absolute noise level added to the simulated flow field.
%                   The standard deviation of the noise is this level 
%                   (which is a percentage) of the median speed.
% simRelNoiseLevel: Relative noise level added to the simulated flow field.
%                   The standard deviation of the noise added to each flow vector 
%                   is this level (which is a percentage) of flow vector length.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Miscellaneous parameters:
%    debugMode: Whether debugging mode is 'on' or 'off'. When it is 'on', some additional results
%               will be calculated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract imaging parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = load([projDir filesep 'lastProjSettings.mat']);
physiParam = s.projSettings.physiParam;

actPixelSize  = 67; %Unit, nm.
frameInterval = 10; %Unit, sec.
if iscell(physiParam)
   %We assume the actin is the 1st image channel.
   actPixelSize  = physiParam{1}.pixelSize;
   frameInterval = physiParam{1}.frameInterval;
end

trackMethod       = 'corr';
isFieldBndFixed   = 'yes';
isDataPosFixed    = 'no';
YModulVariation   = 'homo';
corLen            = 5; %It is supposed to be small and no filtering. Unit, pixels.
sCorLen           = 10; %40; % Used to calculate smoothly interpolated field. Unit, pixels.
gridDx            = 10;  %The grid size (in pixels) in the x-direction.
gridDy            = 10;  %The grid size (in pixels) in the y-direction.
edgCorLen         = 10; %For calculating the displacements on the boundary.
dispType          = 'MFAverage'; %Values: 'SFrame', 'MFAverage' or 'MFTrack'.
showSmooth        = 'yes'; %Values: 'yes' or 'no'.
dataSite          = 'original';
dataToUse         = 'interp'; 
bfDisplaySite     = 'dataSite';
bfFwdOpComputed   = 'none'; 
tfFwdOpComputed   = 'none'; 
forceToIdentify   = 'bf'; 
mcfAngle          = pi/4;
adfAngle          = pi/10;
fixMcfProj        = 'yes';
edgBrkDistTF      = 15;
edgDispIntvLen    = 10;
edgCornerSegCut   = 3;
bspOrderTF        = 4;
bfSigma           = 1e5;
tfSigma           = 1e2;
testTimeStep      = 1;
numBSolsPerFile   = 100;
overwriteBSolBF   = 'yes';
bfScale           = 2e5/10;
tfScale           = 1e3/2;
dispScale         = 15/2;
imgChannel        = 1;
relDispImgFrmNo   = 1;
smDispThreshold   = 0.2;
smForceThreshold  = 0.1;
markMcfBnd        = 'yes';
markAdfBnd        = 'yes';
mcfBndColor       = 'w';
adfBndColor       = 'r';
showFlowVec       = 'yes';
showBdfVec        = 'yes';
showMcfVec        = 'yes';
showAdfVec        = 'yes';
showMixZone       = 'yes';
markMixZone       = 'yes';
markCellEdge      = 'yes';
mixfBndColor      = [0.721 0.721 0.721]; %grey.
cellEdgeColor     = 'r';
imgIRange         = [0 1];
maxSPDToShow      = Inf;
minSPDToShow      = -Inf;
maxBDFToShow      = Inf;
minBDFToShow      = -Inf;
maxMCPToShow      = Inf;
minMCPToShow      = -Inf;
bdfColorDispRange = [0 1]; %Full range.
mcfColorDispRange = [0 1]; %Full range.
mcpColorDispRange = [0 1]; %Full range.
adfColorDispRange = [0 1]; %Full range.
spdColorDispRange = [0 1]; %Full range.
showBdfCBar       = 'yes';
showSpdCBar       = 'yes';
showYModCBar      = 'yes';
spdUnit           = 'umPerMin';
forceToSimu       = 'mcfOnly';
simMCFMagFactor   = 2;
simADCMagFactor   = 2;
simAbsNoiseLevel  = [0 0.00 0.0 0.05 0.10 0.10 0.10 0.10 0.1];
simRelNoiseLevel  = [0 0.05 0.1 0.10 0.10 0.15 0.20 0.25 0.3];
debugMode         = 'off';

param.actPixelSize      = actPixelSize;
param.frameInterval     = frameInterval;
param.trackMethod       = trackMethod;
param.isFieldBndFixed   = isFieldBndFixed;
param.isDataPosFixed    = isDataPosFixed;
param.corLen            = corLen;
param.sCorLen           = sCorLen;
param.gridDx            = gridDx;
param.gridDy            = gridDy;
param.edgCorLen         = edgCorLen;
param.dispType          = dispType;
param.showSmooth        = showSmooth;
param.dataSite          = dataSite;
param.dataToUse         = dataToUse;
param.bfDisplaySite     = bfDisplaySite;
param.bfFwdOpComputed   = bfFwdOpComputed;
param.tfFwdOpComputed   = tfFwdOpComputed;
param.forceToIdentify   = forceToIdentify;
param.mcfAngle          = mcfAngle;
param.adfAngle          = adfAngle;
param.fixMcfProj        = fixMcfProj;
param.edgBrkDistTF      = edgBrkDistTF;
param.edgDispIntvLen    = edgDispIntvLen;
param.edgCornerSegCut   = edgCornerSegCut;
param.bspOrderTF        = bspOrderTF;
param.bfSigma           = bfSigma;
param.tfSigma           = tfSigma;
param.testTimeStep      = testTimeStep;
param.numBSolsPerFile   = numBSolsPerFile;
param.overwriteBSolBF   = overwriteBSolBF;
param.bfScale           = bfScale;
param.tfScale           = tfScale;
param.dispScale         = dispScale;
param.imgChannel        = imgChannel;
param.relDispImgFrmNo   = relDispImgFrmNo;
param.smDispThreshold   = smDispThreshold;
param.smForceThreshold  = smForceThreshold;
param.markMcfBnd        = markMcfBnd;
param.markAdfBnd        = markAdfBnd;
param.mcfBndColor       = mcfBndColor;
param.adfBndColor       = adfBndColor;
param.showFlowVec       = showFlowVec;
param.showBdfVec        = showBdfVec;
param.showMcfVec        = showMcfVec;
param.showAdfVec        = showAdfVec;
param.showMixZone       = showMixZone;
param.markMixZone       = markMixZone;
param.markCellEdge      = markCellEdge;
param.mixfBndColor      = mixfBndColor;
param.cellEdgeColor     = cellEdgeColor;
param.imgIRange         = imgIRange;
param.maxSPDToShow      = maxSPDToShow;
param.minSPDToShow      = minSPDToShow;
param.maxBDFToShow      = maxBDFToShow;
param.minBDFToShow      = minBDFToShow;
param.maxMCPToShow      = maxMCPToShow;
param.minMCPToShow      = minMCPToShow;
param.bdfColorDispRange = bdfColorDispRange;
param.mcfColorDispRange = mcfColorDispRange;
param.mcpColorDispRange = mcpColorDispRange;
param.adfColorDispRange = adfColorDispRange;
param.spdColorDispRange = spdColorDispRange;
param.showBdfCBar       = showBdfCBar;
param.showSpdCBar       = showSpdCBar;
param.showYModCBar      = showYModCBar;
param.spdUnit           = spdUnit;
param.forceToSimu       = forceToSimu;
param.simAbsNoiseLevel  = simAbsNoiseLevel;
param.simRelNoiseLevel  = simRelNoiseLevel;
param.simMCFMagFactor   = simMCFMagFactor;
param.simADCMagFactor   = simADCMagFactor;
param.debugMode         = debugMode;
