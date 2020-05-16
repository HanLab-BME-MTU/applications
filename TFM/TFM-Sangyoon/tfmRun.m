function [] = tfmRun(MDPath)
%% set up
% analysisFolder = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2017-06-29//ChoK1_shRNA_WT_Rescue_FACS_5kPa_006';
if isa(MDPath,'MovieData')
    MD=MovieData.load(MDPath.getFullPath); %,'askUser',false,'askUserChannel',false);
else
    MD=MovieData.load(MDPath,'askUser',false); %,'askUserChannel',false);
end
iPack=  MD.getPackageIndex('TFMPackage');
TFMPack = MD.getPackage(iPack);
status = TFMPack.sanityCheck;
%% Traction reconstruction
curProcess=cell(1,6);
for ii=find(~status)
    curProcess{ii} = TFMPack.getProcess(ii);
    if ~isempty(curProcess{ii})
        if ~curProcess{ii}.success_ || ~curProcess{ii}.updated_ || curProcess{ii}.procChanged_
            curProcess{ii}.run
            MD.save
        end
    end
end
