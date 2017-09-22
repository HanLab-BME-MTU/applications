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
	
	% MDCrop.save();

	%% Collect kinROI movies and volume
    buildPoleKT3DMovies(MDOrig,'package',[], ... 
        'trackingKTPackage',GenericPackage(MDOrig.getPackage(333).processes_(3:4)), ...
        'buildSpindleRefPackage',GenericPackage(MDOrig.getPackage(333).processes_(5:6)) ...
        );
    MDOrig.save();
    
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