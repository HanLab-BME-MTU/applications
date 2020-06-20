function [averagingRange] = findBestCorrelatedRefImage(curRefBeadStack,img)
% function [averagingRange] = findBestCorrelatedRefImage(curRefBeadStack) 
% identifies the best z-location by doing correlation between each z-frame
% with img.
% Sangyoon Han, June 2020

[h,w,numZ] = size(curRefBeadStack);
%Should use only central part
startX=round(w/4); endX=round(3*w/4);
startY=round(h/4); endY=round(3*h/4);

curCor = zeros(numZ,1);
xOffsetAll = zeros(numZ,1);
yOffsetAll = zeros(numZ,1);
% Per image
% midPixelsAllFrames = cell(refMD.zSize_,1);
for jj=1:numZ
    curRefImageFrame = curRefBeadStack(:,:,jj);
    tempRef = curRefImageFrame(startY:endY,startX:endX);
    % Get the cross correlation score
    score = normxcorr2(tempRef,img); 
    % Find peak
    [ypeak, xpeak] = find(score==max(score(:)));
    % Account for the padding that normxcorr2 adds.
    yoffSet = ypeak-size(curRefImageFrame,1);
    xoffSet = xpeak-size(curRefImageFrame,2);
    % Record the max xcorr score
    curCor(jj) = max(score(:));
    xOffsetAll(jj) = xoffSet;
    yOffsetAll(jj) = yoffSet;
end
[~,iMaxCor] = max(curCor);
if iMaxCor>1 && iMaxCor<numZ
    averagingRange = iMaxCor-1:iMaxCor+1;
elseif iMaxCor==1
    averagingRange = iMaxCor:iMaxCor+1;
elseif iMaxCor==numZ
    averagingRange = iMaxCor-1:iMaxCor;
end

end

% This is now obsolete
% % I will take 100th to 200th brightest pixels to guess the best
% % focus (usually the very brightest point is from one
% % extraordinary bead) - SH 20171010
% % Get 100th to 200th pixels
% midPixelsAllFrames = cell(refMD.zSize_,1);
% for jj=1:refMD.zSize_
%     curRefImageFrame = curRefBeadStack(:,:,jj);
%     curRefImageFrameSorted = sort(curRefImageFrame(:),'descend');
%     midPixelsAllFrames{jj} = curRefImageFrameSorted(100:300);
% end
% meanMidInten = cellfun(@mean,midPixelsAllFrames);
% % Take top five frames
% [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
% averagingRange = meanMidIntenIDs(1:5);
% % maxProf = reshape(max(max(curRefBeadStack)),[],1);
% % [~,maxIntenFrame]=max(maxProf);
% % minFocusedFrame=max(1,maxIntenFrame-2);
% % maxFocusedFrame=max(refMD.zSize_,maxIntenFrame+2);
