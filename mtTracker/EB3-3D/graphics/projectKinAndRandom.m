function processList=projectKinAndRandom(MD,processSpindleRef,processDetectPoles,processRandom,kinRange,varargin)
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('name','bundle');
ip.addOptional('processSingleProj',[]);
ip.parse(varargin{:});
p=ip.Results;


%%
relativeOutput='AllKinPrint';

%% loading process results (Cell ref and ISO ref)
kinTracks=load(processSpindleRef.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

tmp=load(processRandom.outFilePaths_{1});
randKinTracks=tmp.randKinTracks;
randKinTracksISO=tmp.randKinTracksISO;
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%% load process indenpendant data
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracksISO=kinTrackData.tracksLabRef;
kinTracksISOInliers=kinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

% Load addPoleRef to weed out inliers kin
poleMovieInfo=load(processDetectPoles.outFilePaths_{1}); poleMovieInfo=poleMovieInfo.poleMovieInfo;
[poleRefsISO,P1,P2]=buildSpindleRef(poleMovieInfo,1);

%%
processList=[];
for kinIdx=kinRange
  kinTrack=kinTracksISOInliers(kinIdx);
  randKinTrack=randKinTracksISOInliers(kinIdx);
  refKP1=buildRefsFromTracks(P1,kinTrack);
  processProj=ExternalProcess(MD,'rawProj');
  project1D(MD,[P1,kinTrack],'dynPoligonREF',[refKP1(1).applyBase(P1,[]) refKP1(1).applyBase(kinTrack,[]) refKP1(1).applyBase(randKinTrack,[])],'FoF',refKP1(1), ...
      'name',[p.name '-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','saveSingleProj',true,'processSingleProj',processProj);
  processList=[processList processProj];
end
