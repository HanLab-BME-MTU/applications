function [meanRefImg,meanRefImgPath] = findBestCorrelatedRefImageMD(refMD,img,iChan)
%[meanRefImg] = findBestCorrelatedRefImageMD(refMD,img,iChan) 
% It is a wrapper function for findBestCorrelatedRefImage. 
% It takes end channel for stack.
% input:
%       refMD           MovieData for image stack of ref img
%       img             one image for current bead image for deformed gel
%       iChan           channel index within refMD if there are more than
%                       one channels
% output:
%       meanRefImg      reference image with the same format as a channel
%       meanRefImgPath  a path for ref img
%   Sangyoon Han June 2020
if nargin<3
    iChan = numel(refMD.channels_);
end
curRefBeadChan = refMD.channels_(iChan);
if refMD.zSize_>1
    curRefBeadStack = curRefBeadChan.loadStack(1);
    averagingRange = findBestCorrelatedRefImage(curRefBeadStack,img);
    curRefBeadStackChosen = curRefBeadStack(:,:,averagingRange);
    meanRefImg = mean(curRefBeadStackChosen,3); 
elseif refMD.nFrames_>1
    curRefBeadStack = zeros(refMD.imSize_(1),refMD.imSize_(2),refMD.nFrames_);
    for ii=1:refMD.nFrames_
        curRefBeadStack(:,:,ii) = curRefBeadChan.loadImage(ii);
    end
    averagingRange = findBestCorrelatedRefImage(curRefBeadStack,img);
    curRefBeadStackChosen = curRefBeadStack(:,:,averagingRange);
    meanRefImg = mean(curRefBeadStackChosen,3); 
else
    meanRefImg = curRefBeadChan.loadImage(1);
end

% store it somewhere
curRef = curRefBeadChan.getPath; %[curRefDir(ii).folder filesep curRefDir(ii).name];
         
[path1,fname1] = fileparts(curRef);
meanRefImgPath = [path1 filesep fname1 '_Chan' num2str(iChan) 'REF.tif'];
imwrite(uint16(meanRefImg),meanRefImgPath,'Compression','none')

end

