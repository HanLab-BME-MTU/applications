function overlayProjTracksList(MD,processProjList,kinTracksInliers,kinTracksISO,randKinTracksInliers,randKinTracksISO,PoleTrack,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('mappedTrackField','associatedMTP1');
ip.addParameter('name','mapped');
ip.parse(MD,varargin{:});
p=ip.Results;
%% loading process results (Cell ref and ISO ref)
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EB3TracksISO=tmp.tracksLabRef;

myColormap=uint8( ...
    [[0 0 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);
% Load addPoleRef to weed out inliers kin
for projIdx=1:length(processProjList)
  process=processProjList(projIdx);
  kinTrack=kinTracksISO(projIdx);
  randKinTrack=randKinTracksISO(projIdx);
  refKP1=buildRefsFromTracks(PoleTrack,kinTrack);
  mappedMT=getfield(kinTracksInliers(projIdx),p.mappedTrackField);
  mappedRef=refKP1.applyBase(EB3TracksISO([mappedMT.index]),[]);
%  mappedSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedMT.index]),[]);

  mappedRandMT=getfield(randKinTracksInliers(projIdx),p.mappedTrackField);
  mappedRandRef=refKP1.applyBase(EB3TracksISO([mappedRandMT.index]),[]);
    %    mappedRandSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]);

  colorIndx=[ones(1,length(mappedMT)) 2*ones(1,length(mappedRandMT))];
  overlayProjTracksMovie(process,'tracks',[mappedRef;mappedRandRef;refKP1.applyBase(kinTrack,[]); refKP1.applyBase(randKinTrack,[])] ,'colorIndx',[colorIndx 3 4],'colormap',myColormap,'name',p.name);
end
