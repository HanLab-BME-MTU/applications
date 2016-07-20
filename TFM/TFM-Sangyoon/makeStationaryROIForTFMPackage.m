function []=makeStationaryROIForTFMPackage(MD)
% function []=makeROImaskForTFMPackage(curMD) creates roiMask.tif for TFM
% package accounting for shift in the bead images
% input: curMD      movieData file containign TFMPackage
% Sangyoon Han 2016. 5.9.
iPack=  MD.getPackageIndex('TFMPackage');
if isempty(iPack)
    error('TFM Package was not run!')
end

curSDCProc = MD.getPackage(iPack).getProcess(1);
if ~isempty(curSDCProc)
    params = curSDCProc.funParams_;
    % Get refROI mask with:
    sampleCellImg = (MD.channels_(2).loadImage(1));
    sampleCellImgD=double(sampleCellImg);
    sampleCellImgNorm=(sampleCellImgD-min(sampleCellImgD(:)))/(max(sampleCellImgD(:))-min(sampleCellImgD(:)));
    beadChannel=MD.channels_(1);
    sampleBeadImg = double(MD.channels_(1).loadImage(MD.nFrames_));
    sampleRefImg = double(imread(params.referenceFramePath));
    compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
    compImg(:,:,3)=sampleCellImgNorm;
    compImg(:,:,2)=(sampleBeadImg-min(sampleBeadImg(:)))/(max(sampleBeadImg(:))-min(sampleBeadImg(:)));
    compImg(:,:,1)=(sampleRefImg-min(sampleRefImg(:)))/(max(sampleRefImg(:))-min(sampleRefImg(:)));
    h1=figure;
    imshow(compImg)
    disp(['Draw rectangle for roi mask for ' MD.movieDataPath_ '. (Red is the reference image.)'])
    if ~isempty(MD.roiMaskPath_)
        roiMask = imread(MD.roiMaskPath_);
        bdryMask=bwboundaries(roiMask);
        hold on
        plot(bdryMask{1}(:,2),bdryMask{1}(:,1),':w')
    end
    if isempty(params.cropROI)
        h=imrect;
    else
        h=imrect(gca,params.cropROI);
    end
%     wait(h);
%     roiMask=createMask(h);
    refROI = wait(h);
    params.cropROI = refROI;
    MD.getPackage(iPack).getProcess(1).setPara(params);
    disp('Done!')
    close(h1)
end