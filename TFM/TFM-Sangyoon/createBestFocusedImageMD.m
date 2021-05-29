function [meanRefImg,meanRefImgPath] = createBestFocusedImageMD(refMD,iChan)
%[meanRefImg] = createBestFocusedImageMD(refMD) It is a wrapper function
%for findBestFocusFromStack. It takes end channel for stack.
%   Sangyoon Han Sep 2018
if nargin<2
    iChan = numel(refMD.channels_);
end
curRefBeadChan = refMD.channels_(iChan);
thresVariance=0.8;
applySobel = true;
if refMD.zSize_>1
    % find the best focus
    curRefBeadStack = curRefBeadChan.loadStack(1);
    averagingRange = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel);
    if numel(averagingRange)>5
        applySobel=false;
        averagingRange = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel);
        if numel(averagingRange)>5
            averagingRange = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel,'amp');
        end
    end
    curRefBeadStackChosen = curRefBeadStack(:,:,averagingRange);
    meanRefImg = mean(curRefBeadStackChosen,3); 
else
    meanRefImg = curRefBeadChan.loadImage(1);
end

% store it somewhere
curRef = curRefBeadChan.getPath; %[curRefDir(ii).folder filesep curRefDir(ii).name];
         
[path1,fname1] = fileparts(curRef);
meanRefImgPath = [path1 filesep fname1 '_Chan' num2str(iChan) '.tif'];
imwrite(uint16(meanRefImg),meanRefImgPath,'Compression','none')

end

