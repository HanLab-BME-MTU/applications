function [stats]=doStackStatistics(stats,pos,doCrop)
if nargin<1
    stats=[];
end

if nargin<2
    pos=length(stats)+1;
end

if nargin<3
    doCrop=0;
end

% if nargin <2 || isempty(inputFileList)
%    [filename, pathname] = uigetfile({'*.TIF;*.tif;*.jpg;*.png;*.*'}, ...
%        'Select First Image to be registered');
%    
%    if ~ischar(filename) || ~ischar(pathname)
%        return;
%    end
%    
%    inputFileList = getFileStackNames([pathname filesep filename]);
% else
%     isValid = 1;
%     for i = 1:numel(inputFileList)
%         isValid = isValid && exist(inputFileList{i}, 'file');
%     end
%     if ~isValid
%         error('Invalid input files.');
%     end
% end

folderNames={'MLC','Actin','Nuclei'};

inputFileListMLC    = getFileListFromFolder('MLC');
inputFileListActin  = getFileListFromFolder('Actin');
inputFileListNuclei = getFileListFromFolder('Nuclei');

for iz = 1:length(inputFileListMLC)
    % Read image.
    Img = double(imread(inputFileListMLC{iz}));
    if iz==1
        ImgSum_MLC=Img;
    else
        ImgSum_MLC=ImgSum_MLC+Img;
    end
end
if doCrop
    [area]=createCropArea(ImgSum_MLC);
else
    area=[1 1 size(ImgSum_MLC)];
end

% imwrite(uint16(R),[target_dir, filesep, 'registered_', name, no,'.tif'],'tiff','Compression','none');

ImgSum_MLC_crop=ImgSum_MLC(area(1):area(3),area(2):area(4));
numPixInROI=numel(ImgSum_MLC_crop);

imagesc(ImgSum_MLC_crop);
axis equal


for iz = 1:length(inputFileListActin)
    % Read image.
    Img = double(imread(inputFileListActin{iz}));
    if iz==1
        ImgSum_Actin=Img;
    else
        ImgSum_Actin=ImgSum_Actin+Img;
    end
end
ImgSum_Actin_crop=ImgSum_Actin(area(1):area(3),area(2):area(4));

for iz = 1:length(inputFileListNuclei)
    % Read image.
    Img = double(imread(inputFileListNuclei{iz}));
    if iz==1
        ImgSum_Nuclei=Img;
    else
        ImgSum_Nuclei=ImgSum_Nuclei+Img;
    end
end
ImgSum_Nuclei_crop=ImgSum_Nuclei(area(1):area(3),area(2):area(4));

ptSpreadFunc=10; %pix

stats(pos).Ippix_MLC_full    = mean(ImgSum_MLC(:));
stats(pos).Ippix_Actin_full  = mean(ImgSum_Actin(:));
stats(pos).Ippix_Nuclei_full = mean(ImgSum_Nuclei(:));

stats(pos).Ippix_MLC    = mean(ImgSum_MLC_crop(:));
stats(pos).Ippix_Actin  = mean(ImgSum_Actin_crop(:));
stats(pos).Ippix_Nuclei = mean(ImgSum_Nuclei_crop(:));

stats(pos).Ippix_MLC_std    = std(   ImgSum_MLC_crop(:));
stats(pos).Ippix_Actin_std  = std( ImgSum_Actin_crop(:));
stats(pos).Ippix_Nuclei_std = std(ImgSum_Nuclei_crop(:));

stats(pos).Ippix_MLC_stdErr95    = facSEMtoSEM95*std(   ImgSum_MLC_crop(:))/sqrt(numPixInROI/ptSpreadFunc.^2);
stats(pos).Ippix_Actin_stdErr95  = facSEMtoSEM95*std( ImgSum_Actin_crop(:))/sqrt(numPixInROI/ptSpreadFunc.^2);
stats(pos).Ippix_Nuclei_stdErr95 = facSEMtoSEM95*std(ImgSum_Nuclei_crop(:))/sqrt(numPixInROI/ptSpreadFunc.^2);

posVec=[100 300 500 700 900];

% 3kPa:
%numNuc=[50 83 108];     % 2011_02_24_IF_A1_pMLC488_actin568_Ecad647_4015_2
%numNuc=[52 72  89 103]; % 2011_02_27_IF_A1_pMLC488_actin568_Ecad647_4015_1
%numNuc=[80 86  96];     % 2011_03_01_IF_A1_pMLC488_actin568_Ecad647_4015

% 65kPa:
%numNuc=[41 91 112 97];    % 2011_02_22_IF_A1_pMLC488_Actin568_Ecad647_1003
%numNuc=[29 89  86 81];    % 2011_03_01_IF_A1_pMLC488_actin568_Ecad647_1003_1
%numNuc=[24 45  72 94 85]; % 2011_03_01_IF_A1_pMLC488_actin568_Ecad647_1003_2

posVec=posVec(1:length(stats));
numNuc=numNuc(1:length(stats));

allStats.stats =stats;
allStats.posVec=posVec;
allStats.numNuc=numNuc;

% Absolute intensities normalized:

figure()
errorbar(posVec,horzcat(stats.Ippix_MLC),horzcat(stats.Ippix_MLC_std));
title('average MLC-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/MLC.fig');

figure()
errorbar(posVec,horzcat(stats.Ippix_Actin),horzcat(stats.Ippix_Actin_std));
title('average Actin-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/Actin.fig');

figure()
errorbar(posVec,horzcat(stats.Ippix_Nuclei),horzcat(stats.Ippix_Nuclei_std));
title('average Nuclei-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/Nuclei.fig');

% MLC-Intensity normalized:

figure()
% the std is:
std_MLC_over_Actin=horzcat(stats.Ippix_MLC)./horzcat(stats.Ippix_Actin).*sqrt((horzcat(stats.Ippix_MLC_std)/horzcat(stats.Ippix_MLC)).^2+(horzcat(stats.Ippix_Actin_std)./horzcat(stats.Ippix_Actin)).^2);
errorbar(posVec,horzcat(stats.Ippix_MLC)./horzcat(stats.Ippix_Actin),std_MLC_over_Actin);
title('tot. MLC-Int. normalized by tot Actin-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/MLCOverActin.fig');

figure()
%std_MLC_over_Nuclei=horzcat(stats.Ippix_MLC)./horzcat(stats.Ippix_Nuclei).*sqrt((horzcat(stats.Ippix_MLC_std)/horzcat(stats.Ippix_MLC)).^2+(horzcat(stats.Ippix_Nuclei_std)./horzcat(stats.Ippix_Nuclei)).^2);
std_MLC_over_Nuclei=horzcat(stats.Ippix_MLC)./horzcat(stats.Ippix_Nuclei).*sqrt((horzcat(stats.Ippix_MLC_std)/horzcat(stats.Ippix_MLC)).^2);
errorbar(posVec,horzcat(stats.Ippix_MLC)./horzcat(stats.Ippix_Nuclei),std_MLC_over_Nuclei);
title('tot. MLC-Int. normalized by tot Nuclei-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/MLCOverNuclei.fig');

figure()
%std_Actin_over_Nuclei=horzcat(stats.Ippix_Actin)./horzcat(stats.Ippix_Nuclei).*sqrt((horzcat(stats.Ippix_Actin_std)/horzcat(stats.Ippix_Actin)).^2+(horzcat(stats.Ippix_Nuclei_std)./horzcat(stats.Ippix_Nuclei)).^2);
std_Actin_over_Nuclei=horzcat(stats.Ippix_Actin)./horzcat(stats.Ippix_Nuclei).*sqrt((horzcat(stats.Ippix_Actin_std)/horzcat(stats.Ippix_Actin)).^2);
errorbar(posVec,horzcat(stats.Ippix_Actin)./horzcat(stats.Ippix_Nuclei),std_Actin_over_Nuclei);
title('tot. Actin-Int. normalized by tot Nuclei-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/ActinOverNuclei.fig');


% MLC-Intensity normalized by number of cells:
%return;

figure()
errorbar(posVec,horzcat(stats.Ippix_MLC_full)./numNuc,facSEMtoSEM95*horzcat(stats.Ippix_MLC_std)./numNuc./sqrt(numNuc));
title('MLC-Int normalized by Nuclei Number');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/MLCOverNucleiNumber.fig');

figure()
errorbar(posVec,horzcat(stats.Ippix_Actin_full)./numNuc,facSEMtoSEM95*horzcat(stats.Ippix_Actin_std)./numNuc./sqrt(numNuc));
title('Actin-Int. normalized by Nuclei Number');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'../figures/ActinOverNucleiNumber.fig');

save('../allStats.mat','allStats');




