function debugEB3Track()
%% a Script/function for systematic refinement of KT tracks

%% loading timing based movie table
allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
blurrPoleCheckedMoviesIdx=(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell));
blurrPoleCheckedMovies=allMovieToAnalyse(blurrPoleCheckedMoviesIdx,:);
goodAndOKSNRIdx=ismember(allMovieToAnalyse.EB3SNR,'OK')|ismember(allMovieToAnalyse.EB3SNR,'Good');
blurrPoleCheckedMoviesHighSNR=allMovieToAnalyse(goodAndOKSNRIdx&blurrPoleCheckedMoviesIdx,:);

% Loading a selected cell of interest
MD=MovieData.loadMatFile(blurrPoleCheckedMovies.analPath{ismember(blurrPoleCheckedMovies.Cell,'cell1_12_halfvol2time')&ismember(blurrPoleCheckedMovies.Setup_min_,'1')});
outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/debugCell12EB3/';

%% Debug Crop
if(isempty(MD.findProcessTag('Crop3D_shorter_Amira_comp','safeCall',true)))
	MDCropRepair=crop3D(MD,(MD.getChannel(1).loadStack(1)),'keepFrame',1:5,'name','shorter_Amira_comp');
	MDCropRepair.sanityCheck();
	MD.save();
else
	MDCropRepair=MovieData.loadMatFile(MD.findProcessTag('Crop3D_shorter_Amira_comp').outFilePaths_{1});
end

% MD=MDCropRepair;
% outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/debugCell12EB3-crop/';

%%
buildAndProjectSpindleRef(MD);
pack=MD.searchPackageName('trackEB3','selectIdx','last');
% pack.eraseProcess(6);
trackEB3(MD,'package',pack,'debug',true,'dynROIView',MD.searchPackageName('dynROIView','selectIdx','last'));
pack=MD.searchPackageName('trackEB3','selectIdx','last');
printProcMIPArray({pack.getProcess(6)},[outputFolder 'tracks'],'MIPIndex',4,'MIPSize',800,'maxWidth',2000);

MD.save();
