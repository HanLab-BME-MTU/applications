function processList=projectKinAndRandom(MD,processSpindleRef,processDetectPoles,processRandom,kinRange,varargin)
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('name','bundle');
ip.addOptional('processSingleProj',[]);
ip.addOptional('showRand',false);
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
%randomMax=processRandom.getParameters().randomDist;
randomMax=20;
processList=[];
for kinIdx=kinRange
  kinTrack=kinTracksISOInliers(kinIdx);
  randKinTrack=randKinTracksISOInliers(kinIdx);
  refKP1=buildRefsFromTracks(P1,kinTrack);
  processProj=ExternalProcess(MD,'rawProj');

  %% build the ROI inset
  insetROI=[P1,kinTrack];
  if p.showRand
      insetROI=[P1,randKinTrack];
  end

  %% Describe the Dynamical ROI ( A pyramid that descibe the maximum distance for randomization)
  baseX=refKP1.getTracksFromBaseVector('X').getMultCoord(randomMax);
  baseY=refKP1.getTracksFromBaseVector('Y').getMultCoord(randomMax);
  % build pyramid
  dynROI=[P1  ...
      kinTrack.getAddCoord(baseX) kinTrack.getAddCoord(baseY) ...
      kinTrack.getAddCoord(baseX.getMultCoord(-1))  kinTrack.getAddCoord(baseY.getMultCoord(-1))];
  arrayfun(@(t) refKP1.applyBase(t,'P1K'),dynROI,'unif',0);

  % Project around the pyramid and show inset.
  project1D(MD,insetROI,'dynPoligonREF',[dynROI.P1K],'FoF',refKP1, ...
    'name',[p.name '-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[1 50]);
  processList=[processList processProj];
end
