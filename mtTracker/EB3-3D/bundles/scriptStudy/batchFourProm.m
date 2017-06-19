ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML1min.mat');
%%
for i=1:ML1min.getSize()
    pack=run1DManifoldDetectorPackage(ML1min.getMovie(i));
    MD.setPackage(333,pack)
    MD.save();
end

%%
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML4min.mat');
for i=1:ML4min.getSize()
    pack=run1DManifoldDetectorPackage(ML4min.getMovie(i));
    MD.setPackage(333,pack)
    MD.save();
end