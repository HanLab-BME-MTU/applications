function [] = changeROIMask(MO)
% function changeROIMask(MD) finds ROI mask associated with the movieData
% MD, ask a user to define the rectangular ROI and save the roiMask.tif
% input: 
%   MD: MovieData or MovieList
% Sangyoon Han, Feb 2018

if isa(MO,'MovieData')
    MO = MovieData.load(MO.getFullPath);
    MDAll{1} = MO;
    nMovies = numel(MDAll);
elseif isa(MO,'MovieList')
    ML = MovieList.load(MO.getFullPath);
    MDAll = ML.movies_;
    nMovies = numel(MDAll);
end
%% ROI for cell to reduce memory burden
for ii=1:nMovies
    curMovie = MDAll{ii};
    nChannels = numel(curMovie.channels_);
    nFrames=curMovie.nFrames_;
    iCellChan = min(2,nChannels);
    if isempty(curMovie.roiMaskPath_) || ~exist(curMovie.roiMaskPath_,'file')
        sampleCellImgFirst = double(curMovie.channels_(iCellChan).loadImage(1));
        sampleCellImgNorm1=(sampleCellImgFirst-min(sampleCellImgFirst(:)))/(max(sampleCellImgFirst(:))-min(sampleCellImgFirst(:)));
        sampleCellImgLast = (curMovie.channels_(iCellChan).loadImage(nFrames));
        sampleCellImgDLast=double(sampleCellImgLast);
        sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
        sampleBeadImg = curMovie.channels_(1).loadImage(1);
        sampleBeadImgD=double(sampleBeadImg);
        sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
        compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
        compImg(:,:,3)=sampleBeadImgNorm;
        compImg(:,:,2)=sampleCellImgNorm2;
        compImg(:,:,1)=sampleCellImgNorm1;
        h12=figure; imshow(compImg), hold on
        disp(['Draw rectangle for ROI for ' curMovie.movieDataPath_ '.'])

        h=imrect;
        ROI_rect = wait(h);
        roiMask=createMask(h);

        % Save it as ROI mask associated with MD
        roiPath=[curMovie.outputDirectory_ filesep 'roiMask.tif'];
        imwrite(roiMask,roiPath);
        curMovie.setROIMaskPath(roiPath);
        % maskArray = imread(curMovie.roiMaskPath_);
        curMovie.roiMask=roiMask;
        close(h12)
    else
        sampleCellImgFirst = (curMovie.channels_(iCellChan).loadImage(1));
        sampleCellImgD1=double(sampleCellImgFirst);
        sampleCellImgNorm1=(sampleCellImgD1-min(sampleCellImgD1(:)))/(max(sampleCellImgD1(:))-min(sampleCellImgD1(:)));
        sampleCellImgLast = (curMovie.channels_(iCellChan).loadImage(nFrames));
        sampleCellImgDLast=double(sampleCellImgLast);
        sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
        sampleBeadImg = curMovie.channels_(1).loadImage(1);
        sampleBeadImgD=double(sampleBeadImg);
        sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
        compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
        compImg(:,:,3)=sampleBeadImgNorm;
        compImg(:,:,2)=sampleCellImgNorm2;
        compImg(:,:,1)=sampleCellImgNorm1;
        h12=figure; imshow(compImg), hold on
        disp(['Draw rectangle for ROI for ' curMovie.movieDataPath_ '.'])
        roiMask=imread(curMovie.roiMaskPath_);
        boundROI=bwboundaries(roiMask);
        hold on, 
        try
            plot(boundROI{1}(:,2),boundROI{1}(:,1),'w--')
            h=imrect;
        catch                
            h=imrect;
        end
        ROI_rect = wait(h);
        roiMask=createMask(h);

        copyfile(curMovie.roiMaskPath_,[curMovie.roiMaskPath_(1:end-4) '_original.tif'])
        % Save it as ROI mask associated with curMovie
        roiPath=[curMovie.outputDirectory_ filesep 'roiMask.tif'];
        imwrite(roiMask,roiPath);
        curMovie.setROIMaskPath(roiPath);
        % maskArray = imread(curMovie.roiMaskPath_);
        curMovie.roiMask=roiMask;
        close(h12)
    end
end
