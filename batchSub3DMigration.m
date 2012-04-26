function batchSub3DMigration(n)

projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/CK_and_FH2_data/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/4D gfpMIIB fix and stain/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/3Dfixset2-2/movieList.mat';

ML = MovieList.load(projPath,0);

% p.BatchMode = true;
% p.CurvSampRad = 3;
% analyze3DMovieMaskedIntensities(ML.movies_{n},p);
% 
% %TEMP - for re-running with new post-processing tracking
% segChans = ML.movies_{n}.processes_{1}.funParams_.ChannelIndex;
% 
% runArgs = {'BatchMode',true,'ChannelIndex',1,'ForceRun',[1 1 1 1 -1],'1ChannelIndex',segChans,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
% 
% process3DMigrationMovie(ML.movies_{n},runArgs{:});

%analyze3DMovieMaskedIntensities(ML.movies_{n});


%  cd('/files/.retain-snapshots.d7d-w0d/LCCB/nih/4D gfpMIIB fix and stain/')
% 
% load('movie array ALL MOVIES Orchestra');
% 
runArgs = {'BatchMode',true,'ChannelIndex',1,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};

process3DMigrationMovie(ML.movies_{n},runArgs{:});



% cd('/files/.retain-snapshots.d7d-w0d/LCCB/nih')
% 
% load('movie array ALL MOVIES')
%
% p.BatchMode = true;
% p.Method = 'SurfaceEnhancement';
% p.PostProcess = true;
% p.FixJumps = false;
% p.ChannelIndex = 1;
% p.ThresholdValue = [];
% p.PreFilterSig = 0;
% 
% MAorchestra(n) = segment3DMovie(MAorchestra(n),p);
% 


% cd('/files/.retain-snapshots.d7d-w0d/LCCB/nih')
% 
% load('movie array ALL MOVIES')
%
% p.BatchMode = true;
% p.Method = 'SurfaceEnhancement';
% p.PostProcess = true;
% p.FixJumps = false;
% p.ChannelIndex = 1;
% p.ThresholdValue = [];
% p.PreFilterSig = 0;
% 
% MAorchestra(n) = segment3DMovie(MAorchestra(n),p);
% 
% 
% p.BatchMode = true;
% p.ForceRun = [-1 1 1 1 1 ];
% p.ChannelIndex = 1;
% 
% process3DMigrationMovie(MAorchestra(n),p);