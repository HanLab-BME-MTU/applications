function [] = runTFMPackage(pathToMovieData)
% runTFMPackage.m is a function for TFM analysis when all set-up is
% already done for movieData in pathToMovieData.
% e.g. analysisFolder = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/2015-07-24/Paxillin1';
MD=MovieData.load(pathToMovieData);
iPack=  MD.getPackageIndex('TFMPackage');
TFM = MD.getPackage(iPack);
%% SDC
prevProcChanged = false;
SDCprocess = TFM.getProcess(1);
if ~SDCprocess.success_ || SDCprocess.procChanged_
    SDCprocess.run
    MD.save
    prevProcChanged = true;
end
%% Displacement
displProcess = TFM.getProcess(2);
if ~displProcess.success_  || displProcess.procChanged_ || prevProcChanged
    displProcess.run
    MD.save
    prevProcChanged = true;
end
%% Displacement correction
displCoProcess = TFM.getProcess(3);
if ~displCoProcess.success_  || displCoProcess.procChanged_ || prevProcChanged
    displCoProcess.run
    MD.save
    prevProcChanged = true;
end
%% Traction reconstruction
forceProcess = TFM.getProcess(4);
if ~forceProcess.success_ || forceProcess.procChanged_ || prevProcChanged
    forceProcess.run
    MD.save
end