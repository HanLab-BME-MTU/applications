function [] = tfmRun(MDPath)
%% set up
% analysisFolder = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2017-06-29//ChoK1_shRNA_WT_Rescue_FACS_5kPa_006';
MD=MovieData.load(MDPath);
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
    prevProcChanged = true;
end
%% Strain Energy quantification
seProcess = TFM.getProcess(5);
if ~seProcess.success_ || seProcess.procChanged_ || prevProcChanged
    seProcess.run
    MD.save
end
