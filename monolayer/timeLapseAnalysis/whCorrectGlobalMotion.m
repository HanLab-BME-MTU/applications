% Function name: whCorrectGlobalMotion
% Author: Assaf Zaritsky, May 2015
% Description: corrects global motion in the image (cause by microscope
%              repeatability problem). Uses the background texture to estimate local 
%              background motion. If most of motion ~= 0 than there is a global motion to be corrected 

function [] = whCorrectGlobalMotion(params,dirs)

% move back to original files
if exist([dirs.mfDataOrig filesep '001_mf.mat'],'file')    
    % unix(sprintf('cp -R %s %s',[dirs.mfDataOrig '*.mat'],dirs.mfData));
    copyfile([dirs.mfDataOrig '*.mat'], dirs.mfData);
end
    
% unix(sprintf('cp -R %s %s',[dirs.mfData '*.mat'],dirs.mfDataOrig));
copyfile([dirs.mfData '*.mat'], dirs.mfDataOrig);


correctionsDx = [];
correctionsDy = [];
medianPrecentDx = [];
medianPrecentDy = [];
% transDx = [];
% transDy = [];

% dilMask = ones(30); % 30 pixels dilation
% dilSe = strel(dilMask);
% 
% [optim,metric]=imregconfig('monomodal');

for t = 1 : params.nTime - params.frameJump
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
    roiFname = [dirs.roiData sprintf('%03d',t) '_roi.mat'];    
    
    fprintf(sprintf('correcting motion estimation frame %d\n',t));
    
    load(mfFname); % dxs, dys
    load(roiFname); % ROI
    
    [correctDx, medianPDx] = getCorrection(dxs,ROI);
    %     xsBackground = dxs(~ROI);
    %     xsBackground = xsBackground(~isnan(xsBackground));
    %     nXsBackground = length(xsBackground);
    %     medianXsBackground = median(xsBackground(:));
    %     if medianXsBackground ~= 0 && (sum(xsBackground == medianXsBackground) > 0.6 * nXsBackground)
    %         correctDx = -medianXsBackground;
    %     else
    %         correctDx = 0;
    %     end
    correctionsDx = [correctionsDx correctDx];
    medianPrecentDx = [medianPrecentDx medianPDx];
    
    [correctDy, medianPDy] = getCorrection(dys,ROI);
    %     ysBackground = dys(~ROI);
    %     ysBackground = ysBackground(~isnan(ysBackground));
    %     nYsBackground = length(ysBackground);
    %     medianYsBackground = median(ysBackground(:));
    %     if medianYsBackground ~= 0 && (sum(ysBackground == medianYsBackground) > 0.6 * nYsBackground)
    %         correctDy = -medianYsBackground;
    %     else
    %         correctDy = 0;
    %     end
    correctionsDy = [correctionsDy correctDy];
    medianPrecentDy = [medianPrecentDy medianPDy];
    
% %     % Compare with imregtform -- does not work because of motion of
% debris / cells!!
% %     imgFname0 = [dirs.images sprintf('%03d',t) '.tif'];
% %     imgFname1 = [dirs.images sprintf('%03d',t+params.frameJump) '.tif'];
% %     roiFname0 = [dirs.roiData sprintf('%03d',t) '_roi.mat'];
% %     roiFname1 = [dirs.roiData sprintf('%03d',t+params.frameJump) '_roi.mat'];
% %     I0 = imread(imgFname0);    
% %     I1 = imread(imgFname1);
% %     load(roiFname0); ROI0 = ROI; clear ROI;
% %     %     load(roiFname1); ROI1 = ROI; clear ROI;
% %     ROI0 = imdilate(ROI0,dilSe);
% %     %     ROI1 = imdilate(ROI1,dilSe);
% %     %     stats0 = regionprops(ROI0,'BoundingBox'); x = stats0.BoundingBox(3);
% %     stats0 = regionprops(ROI0,'BoundingBox'); x = stats0.BoundingBox(3);
% %     I0crop = I0(:,x:end);
% %     I1crop = I1(:,x:end);
% %     [tform]=imregtform(I0crop,I1crop,'translation',optim,metric);
% %     transDx = [transDx -tform.T(3,1)];
% %     transDy = [transDy -tform.T(3,2)];
    
    %     % BELOW: just dilate the ROI A LOT and then take a bounding box
    %     at the ~ROI region
    %     %     I0(ROI0) = nan;
    %     %     I1(ROI1) = nan;
    %     %     I0(ROI0) = 0;
    %     %     I1(ROI1) = 0;
    %     
    %     
    
    %     %     [dydx1, dys1, dxs1, scores1] = blockMatching(I0, I1, params.patchSize,params.searchRadiusInPixels,true(size(I0)),round(correctDx),round(correctDy)); % block width, search radius,
    if abs(correctDx) > 0.5 || abs(correctDy) > 0.5
        imgFname0 = [dirs.images sprintf('%03d',t) '.tif'];
        imgFname1 = [dirs.images sprintf('%03d',t+params.frameJump) '.tif'];
        I0 = imread(imgFname0);    
        I1 = imread(imgFname1);
        
        if size(I0,3) > 1
            tmp = I0(:,:,1) - I0(:,:,2);
            assert(sum(tmp(:)) == 0);
            I0 = I0(:,:,1);
            I1 = I1(:,:,1);
        end
        
        [dydx, dys, dxs, scores] = blockMatching(I0, I1, params.patchSize,params.searchRadiusInPixels,true(size(I0)),round(correctDx),round(correctDy)); % block width, search radius,
    end
    
    %     save(mfFname,'dxs', 'dys','scores');
    
    
    %% correction
    dxs = dxs + correctDx;
    dys = dys + correctDy;
    save(mfFname,'dxs','dys','scores');

end

nCorrected = sum(abs(correctionsDx) > 0.5 | abs(correctionsDy) > 0.5);
nCorrectedDx = sum(abs(correctionsDx) > 0.5);
nCorrectedDy = sum(abs(correctionsDy) > 0.5);
precentCorrected = nCorrected/length(correctionsDx);
save([dirs.correctMotion dirs.expname '_correctMotion.mat'],'correctionsDx','correctionsDy','medianPrecentDy','medianPrecentDx','nCorrected','nCorrectedDx','nCorrectedDy','precentCorrected');%,'transDx','transDy'
end

%%
% medianPDx - precent of patches in the median +- 1 range, should be high!
function [correctDx, medianPDx] = getCorrection(dxs,ROI)
xsBackground = dxs(~ROI);
xsBackground = xsBackground(~isnan(xsBackground));
nXsBackground = length(xsBackground);
medianXsBackground = median(xsBackground(:));

sumMedian0 = sum(xsBackground == (medianXsBackground-1));
sumMedian1 = sum(xsBackground == medianXsBackground);
sumMedian2 = sum(xsBackground == (medianXsBackground+1));
allSumMedian = sumMedian0 + sumMedian1 + sumMedian2;

correctDx = 0;
if (allSumMedian > 0.6 * nXsBackground)
    correctDx = -(((sumMedian0 * (medianXsBackground-1)) + (sumMedian1 * medianXsBackground) + (sumMedian2 * (medianXsBackground+1)))/allSumMedian);
end

medianPDx = allSumMedian/nXsBackground;
end