%% load data
ML=MovieList.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\package\trackProcess\analysis\earlyAndLateML.mat');
MD=ML.getMovie(2);
%% Retrieve +Tips and Pole detections and build a package
processEB3=MD.getPackage(333).getProcess(1);
processPole=MD.getPackage(333).getProcess(5);
processROI=MD.getPackage(333).getProcess(6);

pack=GenericPackage({processEB3 processPole processROI});
spindleEnrichmentVsElev(MD,'package',pack);

%% After first computation one can compute on same package
spindleEnrichmentVsElev(MD,'package',MD.getPackage(400));

