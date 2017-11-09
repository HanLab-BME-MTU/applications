function []=makeROImaskForTFMPackage(curMD)
% function []=makeROImaskForTFMPackage(curMD) creates roiMask.tif for TFM
% package accounting for shift in the bead images
% input: curMD      movieData file containign TFMPackage
% Sangyoon Han 2016. 5.9.
iPack=  curMD.getPackageIndex('TFMPackage');
if isempty(iPack)
    error('TFM Package was not run!')
end

curSDCProc = curMD.getPackage(iPack).getProcess(1);
if ~isempty(curSDCProc) %curSDCProc.procChanged_ || ~curSDCProc.success_
    params = curSDCProc.funParams_;
    % Get refROI mask with:
    sampleCellImg = (curMD.channels_(2).loadImage(1));
    sampleCellImgD=double(sampleCellImg);
    sampleCellImgNorm=(sampleCellImgD-min(sampleCellImgD(:)))/(max(sampleCellImgD(:))-min(sampleCellImgD(:)));
%     beadChannel=curMD.channels_(1);
    sampleBeadImg = double(curMD.channels_(1).loadImage(curMD.nFrames_));
    sampleRefImg = double(imread(params.referenceFramePath));
    compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
    compImg(:,:,3)=sampleCellImgNorm;
    compImg(:,:,2)=(sampleBeadImg-min(sampleBeadImg(:)))/(max(sampleBeadImg(:))-min(sampleBeadImg(:)));
    compImg(:,:,1)=(sampleRefImg-min(sampleRefImg(:)))/(max(sampleRefImg(:))-min(sampleRefImg(:)));
    h1=figure;
    imshow(compImg)
    disp(['Draw rectangle for roi mask for ' curMD.movieDataPath_ '. (Red is the reference image.)'])
    if ~isempty(curMD.roiMaskPath_)
        roiMask = imread(curMD.roiMaskPath_);
        bdryMask=bwboundaries(roiMask);
        hold on
        plot(bdryMask{1}(:,2),bdryMask{1}(:,1),':w')
    end
    h=imrect;
    wait(h);
    roiMask=createMask(h);
    imwrite(roiMask,curMD.roiMaskPath_)
    disp('Done!')
    close(h1)
end