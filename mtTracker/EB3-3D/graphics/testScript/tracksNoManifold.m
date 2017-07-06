% Use printMIP output as projection process and display the raw tracks on
% it
% Abort due to large differences in printMIP format (needs fuse printMIP to
% make it compatible). Better to curb project1D that is still lightly used
% and then augment printMIP to reflect the flexibility.

ML16min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/2DImg_for_16min/analysis/movieList.mat');
MD=ML16min.getMovie(3);

myColormap=uint8( ...
    [[0 255 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);
tic;
processDetectPoles=MD.getPackage(333).getProcess(5);
tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);

processProj=ExternalProcess(MD,'rawProj');
project1D(  MD, ...
            'name','fullMIPNoManifold','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',20,'fringeWidth',60, ...
            'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);        
        
overlayProjTracksMovie(processProj,'tracks',[P1 P2], ... 
            'colorIndx',[1 2],'colormap',myColormap,'name','tracks');
toc;




%%
labFoF=FrameOfRef().setOriginFromTrack(origTracks).genCanonicalBase();
tic;
project1D(  MD,[P1,P2],'FoF',refs(1,2),'dynPoligonREF',refs(1,2).applyBase([P1,P2],''), ...
            'name','CroppedSpindleRef','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',60, ...
            'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99       ]);        
toc;
%%
tic;
project1D(  MD,[P1,P2],'crop','full', ...
            'name','fullSpindleTestFullCrop','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',60, ...
            'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);        
toc;
%%
tic;
project1D(  MD,[P1,P2],'crop','manifold', ...
            'name','fullSpindleTestFullCrop','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',20,'fringeWidth',60, ...
            'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);        
toc;

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
project1D(  MD,[P1,P2],'FoF',refs(1,2),'dynPoligonREF',refs(1,2).applyBase(dynFullPolygon,''), ...
            'name','fullSpindleRef','channelRender','grayRed','saveSingleProj',true, ...
            'processSingleProj',processProjSpindleRef,'intMinPrctil',[1 10],'intMaxPrctil',[100 100]);