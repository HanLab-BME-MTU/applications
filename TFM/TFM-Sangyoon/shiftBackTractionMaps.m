function [] = shiftBackTractionMaps(MD)
% function [] = shiftBackTractionMaps(MD) shift traction maps created
% by TFM package, to the coordinate in deformed state of images.
% input:      MD:    movieData file
% You have to have your TFM package run. Also, you have to have
% transformation file obtained using Alignment/Registration Transform
% Process (transformCreationGUI).
% output: the function will create image sequence of TFM images that have
% the same coordinate as the camera used for CFP channel.

%% Load TFMPackage
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
%% Load SDC process
iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
iBeadChan=1;
if ~isempty(iSDCProc)
    SDCProc=MD.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
end
%% Load traction maps
forceFieldStruc=load(forceFieldProc.outFilePaths_{1});
forceField = forceFieldStruc.forceField;
forceFieldShifted = forceFieldStruc.forceFieldShifted;

tractionMaps=load(forceFieldProc.outFilePaths_{2});
tMap = tractionMaps.tMap;
tMapX = tractionMaps.tMapX;
tMapY = tractionMaps.tMapY;
% band=0;
%% output: replace the existing and back it up
forceFieldPath=forceFieldProc.outFilePaths_{1};
pathFF=fileparts(forceFieldPath);
backupFolder=[pathFF 'Backup'];
tractionImgFolder=[pathFF filesep 'tractionImg10x'];
if exist(backupFolder,'dir')
    ii = 1;
    while exist(backupFolder,'dir')
        backupFolder = [pathFF ' Backup' num2str(ii)];
        ii=ii+1;
    end
end
disp('Backing up the existing files...')
tic; mkdir(backupFolder);
copyfile(pathFF, backupFolder), toc
if ~exist(tractionImgFolder,'dir')
    mkdir(tractionImgFolder)
end
iiformat = ['%.' '3' 'd'];

%% Go one by one
progressText(0,'Re-shifting force map')
for ii=1:nFrames
    % Apply back shift to bead image coordinate
    tform2ref = affine2d([1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
    invTForm2bead = invert(tform2ref);
%     figure, imshow(tMap{ii},[])
    I2 = imwarp(tMap{ii}, invTForm2bead);
    I2x = imwarp(tMapX{ii}, invTForm2bead);
    I2y = imwarp(tMapY{ii}, invTForm2bead);
%     figure, imshow(I2,[])
    % look at curr bead image
    currImage = double(SDCProc.loadChannelOutput(iBeadChan,ii));
    if size(I2,1)<size(currImage,1)
        rowMargin = size(currImage,1)-size(I2,1);
        colMargin = size(currImage,2)-size(I2,2);
        
        I3 = padarray(I2,[rowMargin colMargin],'post');
        I3x = padarray(I2x,[rowMargin colMargin],'post');
        I3y = padarray(I2y,[rowMargin colMargin],'post');
    else
        I3 = I2;
        I3x = I2x;
        I3y = I2y;
    end
    % Crop according to where original image was shifted
    currImageOrg = double(MD.channels_(iBeadChan).loadImage(ii));
    maskROI = false(size(I3));
    ymin=round(maxY+T(ii,1)*(T(ii,1)>0)+1);
    ymax=round(maxY+size(currImageOrg,1)+T(ii,1)*(T(ii,1)<0));
    xmin=round(maxX+T(ii,2)*(T(ii,2)>0)+1);
    xmax=round(maxX+size(currImageOrg,2)+T(ii,2)*(T(ii,2)<0));
    maskROI(ymin:ymax,xmin:xmax)=true;
    forceField(ii) = filterDisplacementField(forceField(ii), maskROI);
    forceFieldShifted(ii) = filterDisplacementField(forceFieldShifted(ii), maskROI);
    xminRef=round(maxX+T(ii,2)+1);
    yminRef=round(maxY+T(ii,1)+1);
    forceField(ii).pos(:,1)=forceField(ii).pos(:,1)-xminRef;
    forceField(ii).pos(:,2)=forceField(ii).pos(:,2)-yminRef;
    forceFieldShifted(ii).pos(:,1)=forceFieldShifted(ii).pos(:,1)-xminRef;
    forceFieldShifted(ii).pos(:,2)=forceFieldShifted(ii).pos(:,2)-yminRef;
    I3(~maskROI)=NaN;
    currIm=imcrop(I3,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImX=imcrop(I3x,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImY=imcrop(I3y,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
%     figure, imshow(currIm,[])
%     currImageCell = double(MD.channels_(2).loadImage(ii));
%     figure, imshow(currImageCell,[])
%     figure, imshow(currIm,[]); hold on, quiver(forceField2(1).pos(1:end,1),forceField2(1).pos(1:end,2),forceField2(1).vec(1:end,1),forceField2(1).vec(1:end,2),'g')
    % I will magnify currIm by ten times (e.g. 125.3 Pa -> 1253)
    tMap{ii}=currIm;
    tMapX{ii}=currImX;
    tMapY{ii}=currImY;
    
    imName = [tractionImgFolder filesep 'TFMap10x' num2str(ii,iiformat) '.tif'];
    imwrite(uint16(currIm*10),imName,'Compression','none');            
    progressText(ii/nFrames,'Re-shifting force map')
end
save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY')
save(forceFieldProc.outFilePaths_{1},'forceField','forceFieldShifted')

disp('Done!')

