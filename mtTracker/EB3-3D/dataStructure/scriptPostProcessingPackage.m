
% MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat');
% MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\laterCell1_12\movieData.mat');
MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');


processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
processSpindleRef=ExternalProcess(MD,'addSpindleRef',@(p) addSpindleRef(p.getOwner(),'processDetectPoles',processDetectPoles,'process',p));
processRandomKinRef=ExternalProcess(MD,'randomizeKinSpindleRef',@(p) randomizeKinSpindleRef(p.getOwner(),'processDetectPoles',processDetectPoles,'process',p));

package=GenericPackage({processDetectPoles processSpindleRef processRandomKinRef});

for pIdx=1:length(package.processes_)
    package.getProcess(pIdx).run();
end
MD.deletePackage(1)
MD.addPackage(package)
%%
package.showGUI()
