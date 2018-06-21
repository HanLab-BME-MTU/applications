function [] = faPackageRun(MDPath)
%% set up
% analysisFolder = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2017-06-29//ChoK1_shRNA_WT_Rescue_FACS_5kPa_006';
if isa(MDPath,'MovieData')
    MD=MovieData.load(MDPath.getFullPath);
else
    MD=MovieData.load(MDPath);
end
iPack=  MD.getPackageIndex('FocalAdhesionPackage');
FAPackage = MD.getPackage(iPack);
prevProcChanged = false;
%% Traction reconstruction
curProcess=cell(1,11);
for ii=1:11
    curProcess{ii} = FAPackage.getProcess(ii);
    if ~isempty(curProcess{ii})
        if ~curProcess{ii}.success_ || curProcess{ii}.procChanged_ || prevProcChanged
            curProcess{ii}.run
            MD.save
            prevProcChanged = true;
        end
    end
end


