
 MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat');
% MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\laterCell1_12\movieData.mat');
% MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');

%%
processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
processDetectPoles.run();
%%
processBuildRef=ExternalProcess(MD,'buildRefs',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'process',p));
processBuildRef.run();
%%
% Simulating FoF selection in the GUI
tmp=load(processBuildRef.outFilePaths_{1});
ref=tmp.refs(1,2); %(Pole1-Pole2)
% Launching projection
processProj=ExternalProcess(MD,'projection',@(p) projectFrameOfRef(MD,ref,'process',p));
processProj.run();

package=GenericPackage({processDetectPoles processBuildRef });
MD.deletePackage(1)
MD.addPackage(package)
%%


