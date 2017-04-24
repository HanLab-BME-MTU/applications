function overlayProjTracksList(MD,processProjList,kinTracksInliers,kinTracksISO,randKinTracksInliers,randKinTracksISO,PoleTrack)
%%

%% loading process results (Cell ref and ISO ref)
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EB3TracksISO=tmp.tracksLabRef;

myColormap=uint8( ... 
    [[0 0 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 0]; ... % kinetochore tracks   
    ]);
% Load addPoleRef to weed out inliers kin
for projIdx=1:length(processProjList)
  process=processProjList(projIdx);
  kinTrack=kinTracksISO(projIdx);
  randKinTrack=randKinTracksISO(projIdx);
  refKP1=buildRefsFromTracks(PoleTrack,kinTrack);

  mappedMT=kinTracksInliers(projIdx).associatedMTP1; 
  mappedRef=refKP1.applyBase(EB3TracksISO([mappedMT.index]),[]);
%  mappedSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedMT.index]),[]);

  mappedRandMT=randKinTracksInliers(projIdx).associatedMTP1;
  mappedRandRef=refKP1.applyBase(EB3TracksISO([mappedRandMT.index]),[]);
    %    mappedRandSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]);

  colorIndx=[ones(1,length(mappedMT)) 2*ones(1,length(mappedRandMT))];
  overlayProjTracksMovie(process,'tracks',[mappedRef;mappedRandRef;refKP1.applyBase(kinTrack,[]); refKP1.applyBase(randKinTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','mapped');
end
