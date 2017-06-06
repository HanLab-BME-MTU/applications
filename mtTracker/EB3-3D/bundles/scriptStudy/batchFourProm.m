ML1min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML1min.mat');
%%
for i=1:ML1min.getSize()
    add1DManifoldDetectorPackage(ML1min.getMovie(i));
end

%%
ML4min=MovieList.load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/analysis/ML4min.mat');
for i=2:ML4min.getSize()
    add1DManifoldDetectorPackage(ML4min.getMovie(i));
end