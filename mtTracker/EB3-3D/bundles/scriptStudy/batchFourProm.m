ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML1min.mat');
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML4min.mat');
MLs=[ML1min ML4min];


%% relocate for testing puposes
ML1min=ML1min.addAnalysisFolder('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/','/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/processPackage/');
ML4min=ML4min.addAnalysisFolder('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/','/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/processPackage/');

%%
ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/processPackage/analysis/ML1min.mat');
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/processPackage/analysis/ML4min.mat');
MLs=[ML1min ML4min];

%%
for mIdx=1:length(MLs)
    for i=1:MLs(mIdx).getSize()
        MD=MLs(mIdx).getMovie(i);
        %pack=run1DManifoldDetectorPackage(MD,'package',MD.getPackage(333));
        tic;
        pack=run1DManifoldDetectorPackage(MD);
        toc;
        MD.setPackage(333,pack)
        MD.save();
    end
end

%% Collect scoring process (process 5) and compare

%%
ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/ML1min.mat');
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/ML4min.mat');
ML8min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/2DImg_for_8min/analysis/movieList.mat');
ML16min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/2DImg_for_16min/analysis/movieList.mat');
MLs=[ML1min ML4min ML8min ML16min];

%%
processCell=cell(1,length(MLs));

for mIdx=1:length(MLs)
    for i=1:2%MLs(mIdx).getSize()
        MD=MLs(mIdx).getMovie(i);
        processCell{mIdx}=[processCell{mIdx} MD.getPackage(333).getProcess(7)];
    end
end
conditionBundleStatistics(processCell,{'1min','4min','8min','16min'},'/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/plot/')


