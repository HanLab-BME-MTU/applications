% External process bears the following features
% - Store process output files and setup figures with minimal coding
% - Store parameter without copying paramater. 
% - Easy transition from function to process sequence
% - Rerun function (experimental)

MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat');

processDetectPoles=ExternalProcess(MD,'detectPoles');
detectPoles(MD,'process',processDetectPoles);

processSpindleRef=ExternalProcess(MD,'addSpindleRef');
addSpindleRef('processDetectPoles',processDetectPoles,'process',processSpindleRef)
 
processTipsBias=ExternalProcess(MD,'MTTipsBias');
MTTipsBias('processSpindleRef',processSpindleRef,'process',processTipsBias,'kinRange',1:5);

processProjection=ExternalProcess(MD,'printAllMTTipsKinPoleRef');
printAllMTTipsKinPoleRef('processMTTipsBias',processTipsBias,'process',processProjection,'kinRange',1:5)

%% Process Storage using classic MD list
MD.reset();
MD.addProcess(processDetectPoles);
MD.addProcess(processSpindleRef);
MD.addProcess(processTipsBias);
MD.addProcess(processProjection);

%% Process Storage using GenericPackage
package=GenericPackage(MD);
package.showGUI();

%%
MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat');
processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p));
processSpindleRef=ExternalProcess(MD,'addSpindleRef',@(p) addSpindleRef(p.getOwner(),'processDetectPoles',processDetectPoles,'process',p));
processTipsBias=ExternalProcess(MD,'MTTipsBias',@(p) MTTipsBias('processSpindleRef',processSpindleRef,'process',p,'kinRange',1:5));
processProjection=ExternalProcess(MD,'printAllMTTipsKinPoleRef',@(p) printAllMTTipsKinPoleRef('processMTTipsBias',processTipsBias,'process',p,'kinRange',1:5));

package=GenericPackage({processDetectPoles processSpindleRef processTipsBias processProjection});
%%
for pIdx=1:length(package.processes_)
    package.getProcess(pIdx).run();
end

%%
package.showGUI()
