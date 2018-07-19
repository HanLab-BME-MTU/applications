function [] = changeROIMask(MD)
% function changeROIMask(MD) finds ROI mask associated with the movieData
% MD, ask a user to define the rectangular ROI and save the roiMask.tif
% input: 
%   MD: MovieData
% Sangyoon Han, Feb 2018

%% ROI for cell to reduce memory burden
nChannels = numel(MD.channels_);
nFrames=MD.nFrames_;
iCellChan = min(2,nChannels);
if isempty(MD.roiMask)
    sampleCellImgFirst = double(MD.channels_(iCellChan).loadImage(1));
    sampleCellImgNorm1=(sampleCellImgFirst-min(sampleCellImgFirst(:)))/(max(sampleCellImgFirst(:))-min(sampleCellImgFirst(:)));
    sampleCellImgLast = (MD.channels_(iCellChan).loadImage(nFrames));
    sampleCellImgDLast=double(sampleCellImgLast);
    sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
    sampleBeadImg = MD.channels_(1).loadImage(1);
    sampleBeadImgD=double(sampleBeadImg);
    sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
    compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
    compImg(:,:,3)=sampleBeadImgNorm;
    compImg(:,:,2)=sampleCellImgNorm2;
    compImg(:,:,1)=sampleCellImgNorm1;
    h12=figure; imshow(compImg), hold on
    disp(['Draw rectangle for ROI for ' MD.movieDataPath_ '.'])

    h=imrect;
    ROI_rect = wait(h);
    roiMask=createMask(h);

    % Save it as ROI mask associated with MD
    roiPath=[MD.outputDirectory_ filesep 'roiMask.tif'];
    imwrite(roiMask,roiPath);
    MD.setROIMaskPath(roiPath);
    % maskArray = imread(MD.roiMaskPath_);
    MD.roiMask=roiMask;
    close(h12)
else
    sampleCellImgFirst = (MD.channels_(iCellChan).loadImage(1));
    sampleCellImgD1=double(sampleCellImgFirst);
    sampleCellImgNorm1=(sampleCellImgD1-min(sampleCellImgD1(:)))/(max(sampleCellImgD1(:))-min(sampleCellImgD1(:)));
    sampleCellImgLast = (MD.channels_(iCellChan).loadImage(nFrames));
    sampleCellImgDLast=double(sampleCellImgLast);
    sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
    sampleBeadImg = MD.channels_(1).loadImage(1);
    sampleBeadImgD=double(sampleBeadImg);
    sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
    compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
    compImg(:,:,3)=sampleBeadImgNorm;
    compImg(:,:,2)=sampleCellImgNorm2;
    compImg(:,:,1)=sampleCellImgNorm1;
    h12=figure; imshow(compImg), hold on
    disp(['Draw rectangle for ROI for ' MD.movieDataPath_ '.'])
    roiMask=MD.roiMask;
    boundROI=bwboundaries(roiMask);
    hold on, 
    try
        plot(boundROI{1}(:,2),boundROI{1}(:,1),'w')
        h=imrect;
    catch                
        h=imrect;
    end
    ROI_rect = wait(h);
    roiMask=createMask(h);

    % Save it as ROI mask associated with MD
    roiPath=[MD.outputDirectory_ filesep 'roiMask.tif'];
    imwrite(roiMask,roiPath);
    MD.setROIMaskPath(roiPath);
    % maskArray = imread(MD.roiMaskPath_);
    MD.roiMask=roiMask;
    close(h12)
end
