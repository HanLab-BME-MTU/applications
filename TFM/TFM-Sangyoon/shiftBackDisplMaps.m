function [] = shiftBackDisplMaps(MD)
% function [] = shiftBackDisplMaps(MD) shift displacement maps created
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
% Load the displField
iDisplFieldProc = 2;
displFieldProc=TFMPackage.processes_{iDisplFieldProc};
iCorDisplFieldProc = 3;
corDisplFieldProc=TFMPackage.processes_{iCorDisplFieldProc};
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
%% Load displ maps
displFieldStruc=load(displFieldProc.outFilePaths_{1});
displField = displFieldStruc.displField;
displFieldShifted = displFieldStruc.displFieldShifted;

corDisplFieldStruc=load(corDisplFieldProc.outFilePaths_{1});
corDisplField = corDisplFieldStruc.displField;
corDisplFieldShifted = corDisplFieldStruc.displFieldShifted;

displMaps=load(displFieldProc.outFilePaths_{2});
dMap = displMaps.dMap;
dMapX = displMaps.dMapX;
dMapY = displMaps.dMapY;

corDisplMaps=load(corDisplFieldProc.outFilePaths_{2});
corDMap = corDisplMaps.dMap;
corDMapX = corDisplMaps.dMapX;
corDMapY = corDisplMaps.dMapY;
% band=0;
%% output: replace the existing and back it up
displFieldPath=displFieldProc.outFilePaths_{1}; corDisplFieldPath=corDisplFieldProc.outFilePaths_{1};
pathDF=fileparts(displFieldPath); pathCorDF=fileparts(corDisplFieldPath);
backupDisplFolder=[pathDF 'Backup']; backupCorDisplFolder=[pathCorDF 'Backup'];
displImgFolder=[pathDF filesep 'displImg10x']; corDisplImgFolder=[pathCorDF filesep 'displImg10x'];

if exist(backupDisplFolder,'dir')
    ii = 1;
    while exist(backupDisplFolder,'dir')
        backupDisplFolder = [pathDF ' Backup' num2str(ii)];
        backupCorDisplFolder = [pathCorDF ' Backup' num2str(ii)];
        ii=ii+1;
    end
end
disp('Backing up the existing files...')
tic; mkdir(backupDisplFolder); mkdir(backupCorDisplFolder);
copyfile(pathDF, backupDisplFolder),copyfile(pathCorDF, backupCorDisplFolder); toc
if ~exist(displImgFolder,'dir')
    mkdir(displImgFolder)
    mkdir(corDisplImgFolder)
end
iiformat = ['%.' '3' 'd'];

%% Go one by one
progressText(0,'Re-shifting displacement map')
for ii=1:nFrames
    % Apply back shift to bead image coordinate
    tform2ref = affine2d([1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
    invTForm2bead = invert(tform2ref);
%     figure, imshow(tMap{ii},[])
    I2 = imwarp(dMap{ii}, invTForm2bead);
    I2x = imwarp(dMapX{ii}, invTForm2bead);
    I2y = imwarp(dMapY{ii}, invTForm2bead);
    I2c = imwarp(corDMap{ii}, invTForm2bead);
    I2xc = imwarp(corDMapX{ii}, invTForm2bead);
    I2yc = imwarp(corDMapY{ii}, invTForm2bead);
%     figure, imshow(I2,[])
    % look at curr bead image
    currImage = double(SDCProc.loadChannelOutput(iBeadChan,ii));
    if size(I2,1)<size(currImage,1)
        rowMargin = size(currImage,1)-size(I2,1);
        colMargin = size(currImage,2)-size(I2,2);
        
        I3 = padarray(I2,[rowMargin colMargin],'post');
        I3x = padarray(I2x,[rowMargin colMargin],'post');
        I3y = padarray(I2y,[rowMargin colMargin],'post');
        I3c = padarray(I2c,[rowMargin colMargin],'post');
        I3xc = padarray(I2xc,[rowMargin colMargin],'post');
        I3yc = padarray(I2yc,[rowMargin colMargin],'post');
    else
        I3 = I2;
        I3x = I2x;
        I3y = I2y;
        I3c = I2c;
        I3xc = I2xc;
        I3yc = I2yc;
    end
    % Crop according to where original image was shifted
    currImageOrg = double(MD.channels_(iBeadChan).loadImage(ii));
    maskROI = false(size(I3));
    ymin=round(maxY+T(ii,1)*(T(ii,1)>0)+1);
    ymax=round(maxY+size(currImageOrg,1)+T(ii,1)*(T(ii,1)<0));
    xmin=round(maxX+T(ii,2)*(T(ii,2)>0)+1);
    xmax=round(maxX+size(currImageOrg,2)+T(ii,2)*(T(ii,2)<0));
    maskROI(ymin:ymax,xmin:xmax)=true;
    displField(ii) = filterDisplacementField(displField(ii), maskROI);
    displFieldShifted(ii) = filterDisplacementField(displFieldShifted(ii), maskROI);
    corDisplField(ii) = filterDisplacementField(corDisplField(ii), maskROI);
    corDisplFieldShifted(ii) = filterDisplacementField(corDisplFieldShifted(ii), maskROI);
    xminRef=round(maxX+T(ii,2)+1);
    yminRef=round(maxY+T(ii,1)+1);
    displField(ii).pos(:,1)=displField(ii).pos(:,1)-xminRef;
    displField(ii).pos(:,2)=displField(ii).pos(:,2)-yminRef;
    displFieldShifted(ii).pos(:,1)=displFieldShifted(ii).pos(:,1)-xminRef;
    displFieldShifted(ii).pos(:,2)=displFieldShifted(ii).pos(:,2)-yminRef;
    corDisplField(ii).pos(:,1)=corDisplField(ii).pos(:,1)-xminRef;
    corDisplField(ii).pos(:,2)=corDisplField(ii).pos(:,2)-yminRef;
    corDisplFieldShifted(ii).pos(:,1)=corDisplFieldShifted(ii).pos(:,1)-xminRef;
    corDisplFieldShifted(ii).pos(:,2)=corDisplFieldShifted(ii).pos(:,2)-yminRef;
    I3(~maskROI)=NaN;
    currIm=imcrop(I3,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImX=imcrop(I3x,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImY=imcrop(I3y,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImc=imcrop(I3c,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImXc=imcrop(I3xc,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    currImYc=imcrop(I3yc,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
%     figure, imshow(currIm,[])
%     currImageCell = double(MD.channels_(2).loadImage(ii));
%     figure, imshow(currImageCell,[])
%     figure, imshow(currIm,[]); hold on, quiver(forceField2(1).pos(1:end,1),forceField2(1).pos(1:end,2),forceField2(1).vec(1:end,1),forceField2(1).vec(1:end,2),'g')
    % I will magnify currIm by ten times (e.g. 125.3 Pa -> 1253)
    dMap{ii}=currIm;
    dMapX{ii}=currImX;
    dMapY{ii}=currImY;

    corDMap{ii}=currImc;
    corDMapX{ii}=currImXc;
    corDMapY{ii}=currImYc;
    
    imName = [displImgFolder filesep 'DMap10x' num2str(ii,iiformat) '.tif'];
    imwrite(uint16(currIm*10),imName,'Compression','none');            
    imNameC = [corDisplImgFolder filesep 'DMap10x' num2str(ii,iiformat) '.tif'];
    imwrite(uint16(currImc*10),imNameC,'Compression','none');            
    progressText(ii/nFrames,'Re-shifting displacement map')
end
save(displFieldProc.outFilePaths_{2},'dMap','dMapX','dMapY')
save(displFieldProc.outFilePaths_{1},'displField','displFieldShifted')

dMap=corDMap;
dMapX=corDMapX;
dMapY=corDMapY;
displField=corDisplField;
displFieldShifted=corDisplFieldShifted;
save(corDisplFieldProc.outFilePaths_{2},'dMap','dMapX','dMapY')
save(corDisplFieldProc.outFilePaths_{1},'displField','displFieldShifted')

disp('Done!')

