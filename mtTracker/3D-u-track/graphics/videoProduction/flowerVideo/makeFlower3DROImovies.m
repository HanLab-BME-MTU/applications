function makeFLower3DROIMovies()
	%% loading timing based movie table
	allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
	blurrPoleCheckedMoviesIdx=(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell));
	blurrPoleCheckedMovies=allMovieToAnalyse(blurrPoleCheckedMoviesIdx,:);
	goodAndOKSNRIdx=ismember(allMovieToAnalyse.EB3SNR,'OK')|ismember(allMovieToAnalyse.EB3SNR,'Good');
	blurrPoleCheckedMoviesHighSNR=allMovieToAnalyse(goodAndOKSNRIdx&blurrPoleCheckedMoviesIdx,:);

	%% Loading a selected cell of interest
	MDOrig=MovieData.loadMatFile(blurrPoleCheckedMovies.analPath{ismember(blurrPoleCheckedMovies.Cell,'cell1_12_halfvol2time')&ismember(blurrPoleCheckedMovies.Setup_min_,'1')});

	%% Debug Crop
	if(isempty(MDOrig.findProcessTag('Crop3D_shorter_Amira_comp','safeCall',true)))
		MDCropRepair=crop3D(MDOrig,(MDOrig.getChannel(1).loadStack(1)),'keepFrame',1:5,'name','shorter_Amira_comp');
		MDCropRepair.sanityCheck();
		MDOrig.save();
	else
		MDCropRepair=MovieData.loadMatFile(MDOrig.findProcessTag('Crop3D_shorter_Amira_comp').outFilePaths_{1});
	end

	%% Estimate kinROI
	%[overlayCell]=fiberTrackabilityAnalysis(MDOrig,'package',GenericPackage(MDOrig.getPackage(333).processes_(1:7)),'forceRunIdx',[],'printManifCount',2,'KT',1);
	MD=MDOrig;
	% MDCrop.save();
	
	buildAndProjectSpindleRef(MD);

	%% Tracks Kinetorchore
	pack=MD.searchPackageName('trackKT','selectIdx','last');
	trackKT(MD,'package',pack,'dynROIView',MD.searchPackageName('dynROIView','selectIdx','last'));
	
	pack=MD.searchPackageName('trackKT','selectIdx','last');
	MD.save();

	%% Build Kin ROI and display only the large lifetime as a sanity check (150)
 	outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/HaysStack/KTDynROI';
    mkdirRobust(outputFolder);
	allKTs=TracksHandle(pack.getProcess(3).loadChannelOutput(2));
	KTs=allKTs([allKTs.lifetime]==MD.nFrames_);
	fillGaps(KTs);
	poleProcess=MD.searchPackageName('buildAndProjectSpindleRef').getProcess(1);
	tmp=load(poleProcess.outFilePaths_{1});
	P1=tmp.tracks(1);
	P2=tmp.tracks(2);	

 	%dynROIs=[build1DDynROI(P1,KTs(1),8) build1DDynROI(P2,KTs(10),8)];
    dynROIs=[build1DDynROI(P1,KTs(1),8,'kinDistCutoff',[-20,10]) build1DDynROI(P2,KTs(10),8,'kinDistCutoff',[-20,10])];

    %% Rendering for ROI selection
% 	dynROIs=renderPKTDynROI(MD,dynROIs);
	% allEBs=TracksHandle(MD.searchPackageName('trackEB3','selectIdx','last').getProcess(3).loadChannelOutput(1));

 %  	tmp=load(fullfile(outputFolder,'dynROIs.mat'));dynROIs=tmp.dynROIs;
 %  	overlayCell=cell(1,numel(dynROIs));
 %  	for KTIndex=1:length(dynROIs)
 %  	    dynROI=dynROIs(KTIndex);
 %  	    processSelectedDetectionOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
 %  	    myColormap=255*flipud(jet(256));
 %  	    [mappedTracks]=mapTracksTo1DManifold(dynROI.ROI,allEBs,dynROI.mappingDistance,'position','endma');
 %  	    overlayProjTracksMovie(dynROI.renderer,'detections', dynROI.ref.applyBase(mappedTracks,''), ... 
 %  	    	'colorIndx',100*ones(size(mappedTracks)), ...
 %  	    	'process',processSelectedDetectionOverlay, ...
 %  	    	'colormap',myColormap,'name','test')
 %  	    overlayCell{KTIndex}=processSelectedDetectionOverlay;
 %  	end;  	
 %    printProcMIPArray(overlayCell,[outputFolder filesep 'tracks'],'MIPIndex',4,'forceSize',false,'MIPSize',400,'maxHeight',1200,'maxWidth',1700);



%     processBuildRef=MD.searchPackageName('buildAndProjectSpindleRef').getProcess(2);
%     dynROIsSpindleRef=build1DDynROI(poleProcess,poleProcess,80);
%     refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
%     dynROIsInSpindleRef=renderPKTDynROI(MD,dynROIs,'contextDynROI',dynROIsSpindleRef(2),'namePrefix','spindleRef','ref',refs(1,2));
%     save(fullfile(outputFolder,'dynROIsInSpindleRef.mat'),'dynROIsInSpindleRef');
%     
%     printProcMIPArray({dynROIsInSpindleRef(1).renderer,dynROIsInSpindleRef(2).renderer,dynROIs(1).renderer,dynROIs(2).renderer},[outputFolder filesep 'candidateKT'],'MIPIndex',4,'forceSize',false,'MIPSize',400,'maxHeight',1200,'maxWidth',1700);

	%% Collect kinROI movies and volume
    % buildPoleKT3DMovies(MD,'package',[], 'name','VisualROI', ... 
    %     'trackingKTPackage',pack, 'dynROIs',dynROIs, ...
    %     'buildSpindleRefPackage',MD.searchPackageName('buildAndProjectSpindleRef') ...
    %     );
    MD.save();
    
    disp('creating subManifold')
    manifold=[P1,KTs(1)];
    kinDistCutoff=[-20,10];
    manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
    manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
    subManifold=[manifold(2).getAddCoord(manifVector.getMultCoord(kinDistCutoff(1)./manifNorm)) manifold(2).getAddCoord(manifVector.getMultCoord(kinDistCutoff(2)./manifNorm))];
    dynROIs=build1DDynROI(subManifold(1),subManifold(2),8,'kinDistCutoff',[-20,10]);
    manifold=[P2,KTs(10)];
    kinDistCutoff=[-20,10];
    manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
    manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
    subManifold=[manifold(2).getAddCoord(manifVector.getMultCoord(kinDistCutoff(1)./manifNorm)) manifold(2).getAddCoord(manifVector.getMultCoord(kinDistCutoff(2)./manifNorm))];
    dynROIs=[dynROIs build1DDynROI(subManifold(1),subManifold(2),8,'kinDistCutoff',[-20,10])];
    
    
    buildPoleKT3DMovies(MD,'package',[], 'name','TestingROI', ... 
        'trackingKTPackage',pack, 'dynROIs',dynROIs, ...
        'buildSpindleRefPackage',MD.searchPackageName('buildAndProjectSpindleRef') ...
        );

    dynROIs=[build1DDynROI(KTs(1),KTs(1),20,'kinDistCutoff',[-20,10]) build1DDynROI(KTs(1),KTs(1),20,'kinDistCutoff',[-20,10])];
 	%% Collect kinROI movies and volume
    buildPoleKT3DMovies(MD,'package',[], 'name','MCMCROI', ... 
        'trackingKTPackage',pack, 'dynROIs',dynROIs, ...
        'buildSpindleRefPackage',MD.searchPackageName('buildAndProjectSpindleRef') ...
        );   
	% %% Seek ROI of interest in Cell of interest
	% pack=MDOrig.getPackage(1001);
	% pack.setProcess(10,[]).setProcess(11,[]).setProcess(12,[]);
	% [overlayCell]=fiberTrackabilityAnalysis(MDOrig,'package',pack,'printManifCount',2,'KT',1:10);
	% MDOrig.save();

	% printProcMIPArray(num2cell([overlayCell{:}]),[outputFolder filesep 'candidateMIPMore'],'MIPIndex',4,'forceSize',false,'MIPSize',200,'maxHeight',1200,'maxWidth',1700);

	% %% Debug trackability (see batchProcessKinROI for cell selection)
	% outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/trackability/1min_cell1_12_halfvol2time/KTDynROI';
	% printProcMIPArray(num2cell([overlayCell{:}]),[outputFolder filesep 'candidateMIPMore'],'MIPIndex',4,'forceSize',false,'MIPSize',200,'maxHeight',1200,'maxWidth',1700);

	% %% print first trackability success
	% outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/trackability/firstSucess/KTDynROI';
	% filePath='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/firstTrackabilitySuccessCell.xlsx';
	% firstTrackabilitySuccessCell=readtable(filePath);
	% for mIdx=1:height(firstTrackabilitySuccessCell)
	%     MD=MovieData.loadMatFile(firstTrackabilitySuccessCell.analPath{mIdx});
	%     [overlayCell]=fiberTrackabilityAnalysis(MD,'package',GenericPackage(MD.getPackage(333).processes_(1:7)),'printManifCount',2,'KT',1:10);
	% end
	% %%
	% printProcMIPArray(num2cell([overlayCell{:}]),[outputFolder filesep 'candidateMIPMore'],'MIPIndex',4,'forceSize',false,'MIPSize',200,'maxHeight',1200,'maxWidth',1700);