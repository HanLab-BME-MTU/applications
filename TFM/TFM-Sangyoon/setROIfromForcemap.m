function [] = setROIfromForcemap(MD)
% This function shows you (FTTC-based) TFM map and lets you draw ROI for
% your BEM (L1 or L2) force reconstruction. FTTC should've been processed
% before selecting this function:

% Load TFM map
tfmPackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load force process
forceProc = tfmPackage.getProcess(4);

% Get force map
tMap = load(forceProc.outFilePaths_{2});
try
    tMap = tMap.tMap;
catch
    % If there is nothing, run the package with FTTC setting
    funParams = forceProc.funParams_;
    funParams.method='FTTC';
    forceProc.setPara(funParams);
    forceProc.run
    tMap = load(forceProc.outFilePaths_{2});
    tMap = tMap.tMap;
end

% Show the map
h2=figure; imshow(tMap{1},[]), colormap jet;

% Draw the ROI
disp(['Draw rectangle for ROI for ' MD.movieDataPath_ '.'])
if isempty(MD.roiMask)
    h=imrect;
else
    h=imrect(gca,MD.roiMask);
end
ROI_rect = wait(h);
roiMask=createMask(h);

% Save it as ROI mask associated with MD
roiPath=[MD.outputDirectory_ filesep 'roiMask.tif'];
imwrite(roiMask,roiPath);
MD.setROIMaskPath(roiPath);
% maskArray = imread(MD.roiMaskPath_);
MD.roiMask=roiMask;
% maskArray = MD.getROIMask;
close(h2)
disp('ROI created!')
















