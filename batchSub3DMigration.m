function batchSub3DMigration(n)

%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Low_mag_data/test_set_5_10_2012/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/CK_and_FH2_data/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/4D gfpMIIB fix and stain/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/4D gfpMIIB fix and stain/movieListControlWholeAndROIs.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/3Dfixset2-2/movieList.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Low_mag_data/4D Low mag1/movieList.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Low_mag_data/New data set3_10min 09_2012/movieList.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Low_mag_data/movieListROIsandUncropped.mat';
projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Hunter data 2012_09_fixedcells/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Hunter data 2012_09_fixedcells/phall488-mem/movieListAll488.mat';
%projPath =  '/files/.retain-snapshots.d7d-w0d/LCCB/nih/MyoII-GFP timelapse/4D MyoIIA/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/MyoII-GFP timelapse/4D MyoIIA/movieListROIs.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/lowresIIA/movieListTimeLapseOnly.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/MyoII-GFP timelapse/4D MyoIIA/movieListPhotoBleachCorrROIs.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/myoIIA-20X_60X_set3_20130122/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/myoIIA-20X_60X_set3_20130122/60x/movieListROIs.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/act-mem_2013_01/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/myoIIA60X_2103_02/movieListAllROIs.mat';
%projPath = '/files/.retain-snapshots.%d7d-w0d/LCCB/nih/myoIIA60X_2103_02/movieListAll.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Spiders/movieList.mat';
%projPath = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/movieListSuccessfullyProcessed.mat';


ML = MovieList.load(projPath,0);

% p.BatchMode = true;
% p.CurvSampRad = 3;
% analyze3DMovieMaskedIntensities(ML.movies_{n},p);
% 
% %TEMP - for re-running with new post-processing tracking
% segChans = ML.movies_{n}.processes_{1}.funParams_.ChannelIndex;
% 

%segChans = 2; %Use 560 (tdtomato CAAX) to seg for act+mem data


%runArgs = {'BatchMode',true,'ChannelIndex',1,'ForceRun',[0 1 1 1 1],'1ChannelIndex',segChans,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0,'2SampRad',2e3};
%runArgs = {'BatchMode',true,'ForceRun',[-1 1 1 1 1],'2SampRad',2e3};
%runArgs = {'BatchMode',true,'ChannelIndex',1,'ForceRun',[1 1 1 1 1],'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
runArgs = {'BatchMode',true,'ChannelIndex',1,'1ChannelIndex',1:2,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
%runArgs = {'BatchMode',true,'ChannelIndex',1,'ForceRun',[1 1 1 1 1],'1ChannelIndex',1:2,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
%runArgs = {'BatchMode',true,'ChannelIndex',1,'1ChannelIndex',2,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0,'2SampRad',2e3};


process3DMigrationMovie(ML.movies_{n},runArgs{:});


p.BatchMode = true;
%calc3DMoviePhotobleaching(ML.movies_{n},p)
p.CurvSampRad = 2e3;
%p.PhotoBleachMeth = 'None';
%p.TrendRemoval = 'None';
p.PhotoBleachMeth = 'SelfCortical';
p.TrendRemoval = 'Linear';
%p.ChannelIndex = 1:2;
%p.ChannelIndex = 1; 
%Just use default, which is all channels

analyze3DMovieMaskedIntensities(ML.movies_{n},p);




% load('/files/.retain-snapshots.d7d-w0d/LCCB/nih/movieListSuccesfullyProcessed and indices.mat','isWT','isBleb')
% pBleb.OutputDirectory = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Post Processing Myo Inhib and WT/Mask Geometry/Blebbistatin';
% pWT.OutputDirectory = '/files/.retain-snapshots.d7d-w0d/LCCB/nih/Post Processing Myo Inhib and WT/Mask Geometry/WT';
% if n == 1
%     postProcess3DMovieArrayMaskGeometry([ML.movies_{isBleb}],pBleb);
% else
%     postProcess3DMovieArrayMaskGeometry([ML.movies_{isWT}],pWT);
% end




%  cd('/files/.retain-snapshots.d7d-w0d/LCCB/nih/4D gfpMIIB fix and stain/')
% 
% load('movie array ALL MOVIES Orchestra');
% 

%runArgs = {'BatchMode',true,'ChannelIndex',1,'ForceRun',[1 1 1 1 1],'1ChannelIndex',2,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
%runArgs = {'BatchMode',true,'ChannelIndex',1,'1ChannelIndex',2,'1Method','SurfaceEnhancement','1PostProcess',true,'1FixJumps',false,'1ThresholdValue',[],'1PreFilterSig',0};
%process3DMigrationMovie(ML.movies_{n},runArgs{:});



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