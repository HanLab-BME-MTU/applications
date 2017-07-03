ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML1min.mat'%);
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML4min.mat');

MLs=[ML1min ML4min];

%%
for mIdx=1:length(MLs)
    for i=1:MLs(mIdx).getSize()
        MD=MLs(mIdx).getMovie(i);
        %pack=run1DManifoldDetectorPackage(MD,'package',MD.getPackage(333));
        pack=run1DManifoldDetectorPackage(MD);
        MD.setPackage(333,pack)
        MD.save();
    end
end

%% Collect scoring process (process 5) and compare
processCell=cell(1,2);

for mIdx=1:length(MLs)
    for i=1:MLs(mIdx).getSize()
        MD=MLs(mIdx).getMovie(i);
        processCell{mIdx}=[processCell{mIdx} MD.getPackage(333).getProcess(5)];
    end
end

conditionBundleStatistics(processCell)

