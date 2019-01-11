function []=createTractionMapTiff(MD)
% function []=createTractionMap(MD) creates a heatmap of traction magnitude
% in the coordinates of undeformed position (so that they can be used for
% other applications that does not need 
% without exaggerating it with a constant. It stores the map TIF image to
% /TFMPackage/TractionMapTiff.
% Sangyoon Han 2019 January

%% Load TFMPackage
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));

%% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
pathFF=TFMPackage.outputDirectory_;
tractionImgFolder=[pathFF filesep 'TractionMapTiff'];
if ~exist(tractionImgFolder,'dir')
%     system(['mkdir -p ' tractionImgFolder]);
    mkdir(tractionImgFolder);
end

tMap=forceFieldProc.loadChannelOutput('output','tMapUnshifted'); %'tMap');% in Pa per pixel (1pix x 1pix)
padZeros=floor(log10(nFrames))+1;
%% Image writing
progressText(0,'creating traction mag tiff image') % Update text

for ii=1:nFrames
    %imwrite(uint16(Mblue),[target_dir,filesep,'Force_Magnitude',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none') 
    curTmap = tMap(:,:,ii);
    imwrite(uint16(curTmap),[tractionImgFolder,filesep,'ForceMag',num2str(ii,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none');
    progressText(ii/nFrames,'creating traction mag tiff image') % Update text
end
disp('Tiff TFM map creation done!')
disp(['Saved in ' tractionImgFolder])