function [meanRefImg,meanRefImgPath] = createBestFocusedImageMD(refMD)
%[meanRefImg] = createBestFocusedImageMD(refMD) It is a wrapper function
%for findBestFocusFromStack. It takes end channel for stack.
%   Sangyoon Han Sep 2018
curRefBeadChan = refMD.channels_(end);
thresVariance=0.8;applySobel=true;
if refMD.zSize_>1
    % find the best focus
    curRefBeadStack = curRefBeadChan.loadStack(1);
    averagingRange = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel);
    curRefBeadStackChosen = curRefBeadStack(:,:,averagingRange);
    meanRefImg = mean(curRefBeadStackChosen,3); 
else
    meanRefImg = curRefBeadChan.loadImage(1);
end

% store it somewhere
curRef = curRefBeadChan.getPath; %[curRefDir(ii).folder filesep curRefDir(ii).name];
         
[path1,fname1] = fileparts(curRef);
meanRefImgPath = [path1 filesep fname1 '.tif'];
imwrite(uint16(meanRefImg),meanRefImgPath,'Compression','none')

end

