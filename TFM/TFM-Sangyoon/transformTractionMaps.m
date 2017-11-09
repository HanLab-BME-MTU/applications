function [] = transformTractionMaps(MD,pathToTransformData)
% function [] = transformTractionMaps(MD) transforms traction maps created
% by TFM package, based on transformation acquired by Biosensor package
% input:      MD:    movieData file
%               pathToTransformData: path to transformation file from bead
%               channel to CFP channel.
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
else
    T = zeros(nFrames,2);
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
end
%% Load displacementField
% iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);     
% if isempty(iDisplFieldProc)
%     iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
% end
% displFieldProc=MD.processes_{iDisplFieldProc};
% displField=displFieldProc.loadChannelOutput;
display('Backing up existing traction map...')
tic
[pathstr2,name2,ext2] = fileparts(forceFieldProc.outFilePaths_{2});
[pathstr1,name1,ext1] = fileparts(forceFieldProc.outFilePaths_{1});
backUpPathForceField =  [pathstr1 filesep name1,'BeforeTransform' ext1];
backUpPathTractionMap =  [pathstr2 filesep name2,'BeforeTransform' ext2];
if ~exist(backUpPathTractionMap,'file') && exist(forceFieldProc.outFilePaths_{2},'file')
    copyfile(forceFieldProc.outFilePaths_{1},backUpPathForceField)
    copyfile(forceFieldProc.outFilePaths_{2},backUpPathTractionMap)
end
toc
%% Load traction maps
forceField=forceFieldProc.loadChannelOutput;
tractionMaps=load(forceFieldProc.outFilePaths_{2});
tMap = tractionMaps.tMap;
tMapX = tractionMaps.tMapX;
tMapY = tractionMaps.tMapY;
band=0;
%% Load transform
TrToCFT = load(pathToTransformData);
TrToCFT = TrToCFT.xForm;
n = MD.imSize_(1);
m = MD.imSize_(2);
%% output
calibratedTFMap = [MD.outputDirectory_ filesep 'calibratedTFMap'];
if exist(calibratedTFMap,'dir')
    backupFolder = [calibratedTFMap ' Backup']; % name]);
    ii = 1;
    while exist(backupFolder,'dir')
        backupFolder = [calibratedTFMap ' Backup ' num2str(ii)];
        ii=ii+1;
    end
    mkdir(backupFolder);
    copyfile(calibratedTFMap, backupFolder,'f')
else
    mkdir(calibratedTFMap)
end
iiformat = ['%.' '3' 'd'];
newForceField=forceField;

progressText(0,'Transforming force map')
%% Go one by one
for ii=1:nFrames
    if ~isempty(iSDCProc)
        % Apply back shift to bead image coordinate
        tform2ref = affine2d([1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        invTForm2bead = invert(tform2ref);
    %     figure, imshow(tMap{ii},[])
        I2 = imwarp(tMap{ii}, invTForm2bead);
        I2x = imwarp(tMapX{ii}, invTForm2bead);
        I2y = imwarp(tMapY{ii}, invTForm2bead);
        
        newForceField(ii).pos(:,1) = forceField(ii).pos(:,1)-T(ii, 1);
        newForceField(ii).pos(:,2) = forceField(ii).pos(:,2)-T(ii, 2);
    %     figure, imshow(I2,[])
        % look at curr bead image
        currImage = double(SDCProc.loadChannelOutput(iBeadChan,ii));
        % if the size of I2 is smaller than currImage, pad zeros to the right
        % and bottom
    %     reg_grid1=createRegGridFromDisplField(forceField,1,0); %2=2 times fine interpolation
    %     imSizeX = reg_grid1(end,end,1)-reg_grid1(1,1,1)+1;
    %     imSizeY = reg_grid1(end,end,2)-reg_grid1(1,1,2)+1;
    %     w = imSizeX;
    %     h = imSizeY;
    %     centerX = ((reg_grid1(end,end,1)+reg_grid1(1,1,1))/2);
    %     centerY = ((reg_grid1(end,end,2)+reg_grid1(1,1,2))/2);
    %     xmin = ceil(centerX-w/2+band);
    %     xmax = floor(centerX+w/2-band);
    %     ymin = ceil(centerY-h/2+band);
    %     ymax = floor(centerY+h/2-band);
    %     cropInfo = [xmin,ymin,xmax,ymax];
    %     refFrame = double(imread(SDCProc.outFilePaths_{2,1}));
        % Something weird happened to tMap..
        % There was a bug in calculateMovieForceField that inappropriately
        % refelct SDC-ed reference image. It actually used the size of un-SDC-ed
        % reference image. This bug is removed now (2/16/16), but for nadia's
        % stuff, we need to be able to use the incorrectly adjusted image.
        % ------- THIS PART SHOULD BE USED TEMPORARILLY AND SHOULD BE REMOVED
        % IN THE FUTURE  --------------------------------------- %
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
        I4=imcrop(I3,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
        I4x=imcrop(I3x,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
        I4y=imcrop(I3y,[maxX+T(ii,2)+1,maxY+T(ii,1)+1,size(currImageOrg,2)-1,size(currImageOrg,1)-1]);
    %     figure, imshow(I4,[])
    else
        I4=tMap{ii};
        I4x=tMapX{ii};
        I4y=tMapY{ii};
    end
    %% Apply transform w/ bead image
    currIm = imtransform(I4,TrToCFT,'XData',[1 m],'YData',[1 n],'FillValues',0);
    currImX = imtransform(I4x,TrToCFT,'XData',[1 m],'YData',[1 n],'FillValues',NaN);
    currImY = imtransform(I4y,TrToCFT,'XData',[1 m],'YData',[1 n],'FillValues',NaN);
    tMap{ii} = currIm;
    tMapX{ii} = currImX;
    tMapY{ii} = currImY;
    
    % filter out vectors outside the image
    disp('Filter out forcefield outside the mask...')
    tic
    newForceField(ii)=filterDisplacementField(newForceField(ii),~isnan(currIm));
    toc
    % Re-reconstruct newForceField from currImX and currImY
    for jj=1:length(newForceField(ii).pos)
        newForceField(ii).vec(jj,:)=[currImX(round(newForceField(ii).pos(jj,2)),round(newForceField(ii).pos(jj,1))) ,....
            currImY(round(newForceField(ii).pos(jj,2)),round(newForceField(ii).pos(jj,1)))];
    end
    % I will magnify currIm by ten times (e.g. 125.3 Pa -> 1253)
    imName = [calibratedTFMap filesep 'calibratedTFMap10x' num2str(ii,iiformat) '.tif'];
    imwrite(uint16(currIm*10),imName,'Compression','none');            
    progressText(ii/nFrames,'Transforming force map')

end
forceField=newForceField;
forceFieldShifted=newForceField;
save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY'); 
save(forceFieldProc.outFilePaths_{1},'forceField','forceFieldShifted');

disp('Done!')

