ML=MovieList.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\package\trackProcess\analysis\earlyAndLateML.mat');
MLs=[ML];
MD=ML.getMovie(2);
%% Retrieve +Tips and Pole detections and build a package
processEB3=MD.getPackage(333).getProcess(1);
processPole=MD.getPackage(333).getProcess(5);
processROI=MD.getPackage(333).getProcess(6);

pack=GenericPackage({processEB3 processPole processROI});
spindleEnrichmentVsElev(MD,'package',pack);

%% After first computation one can compute on same package
spindleEnrichmentVsElev(MD,'package',MD.getPackage(400));

%%
processCell=cell(1,length(MLs));

for mIdx=1:length(MLs)
    for i=2%MLs(mIdx).getSize()
        MD=MLs(mIdx).getMovie(i);
        processCell{mIdx}=[processCell{mIdx} MD.getPackage(400).getProcess(5)];
    end
end
conditionEnrichmentStatistics(processCell,{'1min','4min','8min','16min'},'/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/plot/')

