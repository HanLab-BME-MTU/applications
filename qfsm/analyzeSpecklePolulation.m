function [] = analyzeSpecklePolulation(MD)
%unction [] = analyzeSpecklePolulation(MD) analyzes sub-populations of
%stationary vs mobile speckle trajectories by loading them from QFSM
%pakcage of the MD movieData.
%   input:
%           MD  MovieData
%   output is going to be stored ...
% Sangyoon Han and Elizabeth Kaechele
% October 8 2019

%% Load QFSM Package
iQFSM = MD.getPackageIndex('QFSMPackage');
qFSM = MD.getPackage(iQFSM);
%% Get all trajectories
iSpeckleTrackingProcess = 6;
speckleTrackingProcess = qFSM.getProcess(iSpeckleTrackingProcess);

tracks = speckleTrackingProcess.loadChannelOutput(1,'output','MPM');
%         MPM        : Magic Position Matrix 
%         MPM = [ y  x  y  x  y  x ... ]
%                  t1    t2    t3

iFrame=1;
index = all(tracks(:,2*iFrame-1:2*iFrame),2)~=0;
varargout{i}=tracks(index,1:2*iFrame);    

%% Get cell segmentation

%% Separate cell segmentation with LA and LP (via distance from the edge)

%% Using moving window to process the tracks. (related to final project)

%%

end

