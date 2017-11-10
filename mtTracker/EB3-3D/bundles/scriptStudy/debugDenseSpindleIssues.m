% Debug issues with denser MT
% -- projection issues
% -- association issues 
% -- tracking problem

myColormap=uint8( ...
    [[0 255 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);

%% projection issue

ML16min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/2DImg_for_16min/analysis/movieList.mat');
MD=ML16min.getMovie(2);

%% Observe issue
processDetectPoles=MD.getPackage(333).getProcess(5);
processScoring=MD.getPackage(333).getProcess(7);


tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);

tmp=load(processScoring.outFilePaths_{1});
zScores=tmp.zScores;
mappingDist=tmp.mappingDist;  
kinTracksISOInliers=tmp.kinTracksISOInliers;

% issue
tIdx=35;
track=kinTracksISOInliers(tIdx);
refP1P2=buildRefsFromTracks(P1,P2);
ROI=[P1 track];
processProj=ExternalProcess(MD,'projP1');
project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P2 ROI],''),'FoF',refP1P2, ...
    'name',['SP-P1-' num2str(tIdx) '-debug'], ...
    'channelRender','grayRed','processSingleProj',processProj, ...
    'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
overlayProjTracksMovie(processProj,'tracks',[refP1P2.applyBase([track.associatedMTP1; track],[]) ] , ...
    'colorIndx',[ones(size(track.associatedMTP1)); 3],'colormap',myColormap,'name','test');

%% No ref proj








%% Refine project1D to allow for FoF changing at the cells scale
tic;
origTracks=TracksHandle();
origTracks.x=ones(1,MD.nFrames_);
origTracks.y=ones(1,MD.nFrames_);
origTracks.z=ones(1,MD.nFrames_);
origTracks.startFrame=1;
origTracks.endFrame=MD.nFrames_;
maxTracks=TracksHandle();
maxTracks.x=MD.imSize_(2)*ones(1,MD.nFrames_);
maxTracks.y=MD.imSize_(1)*ones(1,MD.nFrames_);
maxTracks.z=ones(1,MD.nFrames_)*MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_;
maxTracks.startFrame=1;
maxTracks.endFrame=MD.nFrames_;
toc;
dynFullPolygon=[origTracks maxTracks];

processDetectPoles=MD.getPackage(333).getProcess(5);
tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);

processBuildRef=MD.getPackage(333).getProcess(6);
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
processProjSpindleRef=ExternalProcess(MD,'rawProj');

%%
project1D(  MD,[P1,P2],'FoF',refs(1,2),'dynPoligonREF',refs(1,2).applyBase(dynFullPolygon,''), ...
            'name','fullSpindleRef','channelRender','grayRed','saveSingleProj',true, ...
            'processSingleProj',processProjSpindleRef,'intMinPrctil',[1 10],'intMaxPrctil',[100 100]);
%%
labFoF=FrameOfRef().setOriginFromTrack(origTracks).genCanonicalBase();
tic;
project1D(  MD,[P1,P2],'FoF',[],'dynPoligonREF',dynFullPolygon, ...
            'name','fullSpindle','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',60, ...
            'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99       ]);        
toc;