function processList=projectKinAndRandom(MD,procSpindleRefOrKin,procDetectPoles,procRandomOrTracks,varargin)
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('kinRange',[]);
ip.addOptional('name','bundle');
ip.addOptional('showRand',false);
ip.parse(varargin{:});
p=ip.Results;

relativeOutput='AllKinPrint';

%% loading process results (Cell ref and ISO ref)
if(isa(procSpindleRefOrKin,'Process'))
  kinTracks=load(procSpindleRefOrKin.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
  kinTracksInliers=load(procSpindleRefOrKin.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

  outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
  kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
  kinTracksISO=kinTrackData.tracksLabRef;
  kinTracksISOInliers=kinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
else
  kinTracksISOInliers=procSpindleRefOrKin;
end

if(isa(procRandomOrTracks,'Process'))
  tmp=load(procRandomOrTracks.outFilePaths_{1});
  randKinTracksISO=tmp.randKinTracks;
  randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
else
  randKinTracksISOInliers=procRandomOrTracks;
end
kinRange=p.kinRange;
if(isempty(kinRange))
  kinRange=1:length(kinTracksISOInliers)
end

% Load addPoleRef to weed out inliers kin
poleMovieInfo=load(procDetectPoles.outFilePaths_{1}); poleMovieInfo=poleMovieInfo.poleMovieInfo;
[poleRefsISO,P1,P2]=buildSpindleRef(poleMovieInfo,1);

%%
%randomMax=processRandom.getParameters().randomDist;
randomMax=20;
processList=[];
for kinIdx=kinRange
  kinTrack=kinTracksISOInliers(kinIdx);
  refKP1=buildRefsFromTracks(P1,kinTrack);
  processProj=ExternalProcess(MD,'rawProj');

  %% build the ROI inset
  insetROI=[P1,kinTrack];
  if p.showRand
      randKinTrack=randKinTracksISOInliers(kinIdx);
      insetROI=[P1,randKinTrack];
  end

  %% Describe the Dynamical ROI ( A pyramid that descibe the maximum distance for randomization)
  baseX=refKP1.getTracksFromBaseVector('X').getMultCoord(randomMax);
  baseY=refKP1.getTracksFromBaseVector('Y').getMultCoord(randomMax);
  % build pyramid
  dynROI=[P1  ...
      kinTrack.getAddCoord(baseX) kinTrack.getAddCoord(baseY) ...
      kinTrack.getAddCoord(baseX.getMultCoord(-1))  kinTrack.getAddCoord(baseY.getMultCoord(-1))];
  dynROIManifRef=arrayfun(@(t) refKP1.applyBase(t,'P1K'),dynROI,'unif',0);

  if(p.showRand)
      name=[p.name '-P1-kin-' num2str(kinTrack.index) '-R'];
  else
      name=[p.name '-P1-kin-' num2str(kinTrack.index)];
  end
  % Project around the pyramid and show inset.
  project1D(MD,insetROI,'dynPoligonREF',[dynROIManifRef{:}],'FoF',refKP1, ...
    'name',name,'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[1 50]);
  processList=[processList processProj];
end
