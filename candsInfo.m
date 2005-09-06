function [candsPerFrame,candsCoords,avg,sigma] = candsInfo
%candsInfo queries user for cands directory and returns statistics
%
%SYNOPSIS   [candsPerFrame,candsCoords,avg,sigma] = candsInfo
%
%INPUT      user is queried (by searchFiles2) for directory where cands 
%           files are located
%
%OUTPUT     candsPerFrame: column vector, nth entry of which is the number
%                          of candidate speckles in file n
%           candsCoords: array containing x,y coordinates for local max and
%                        3 local mins for each candidate speckle: speckles 
%                        X coords (8) X frames
%           avg: mean number of candidate speckles per file, computed over 
%                n files
%           sigma: standard deviation of speckle number over the n frames
%
%c: 9-05 kathryn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% listOfFiles contains only cands***.mat files
[listOfCanstFiles,candsDirectory] = searchFiles2('cands','.mat',[],'ask',0,'all','Select cands directory to evaluate');
numFiles = length(listOfCanstFiles);
str1 = sprintf('Directory searched: %s \r',candsDirectory);
str2 = sprintf('Number of cands files: %d \r',numFiles);
disp(str1);
disp(str2);

avg = 0;
avg = 0;
candsPerFrame = zeros(numFiles,1);

% put number of speckles per cands file in vector candsPerFrame
% ifelse used to accommodate cands file numbering system
for i=1:numFiles
    if i<10
        load([candsDirectory sprintf('\\cands00%d.mat',i)])
        candsPerFrame(i) = length(cands);
    elseif i>=10 & i<100
        load([candsDirectory sprintf('\\cands0%d.mat',i)])
        candsPerFrame(i) = length(cands);
    elseif i>=100
        load([candsDirectory sprintf('\\cands%d.mat',i)])
        candsPerFrame(i) = length(cands);
    else
        disp('function cannot process >= 1000 files');
    end
end

% array containing x,y coordinates for local max and 3 local mins for 
% each candidate speckle: speckles X coords X frames
r = max(candsPerFrame);
candsCoords = -1*ones(r,8,numFiles);
for i=1:numFiles
    if i<10
        load([candsDirectory sprintf('\\cands00%d.mat',i)])
        r = candsPerFrame(i);
        for j=1:r
            candsCoords(j,1:2,i) = cands(j).Lmax;
            candsCoords(j,3:4,i) = cands(j).Bkg1;
            candsCoords(j,5:6,i) = cands(j).Bkg2;
            candsCoords(j,7:8,i) = cands(j).Bkg3;
        end
    elseif i>=10 & i<100
        load([candsDirectory sprintf('\\cands0%d.mat',i)])
        r = candsPerFrame(i);
        for j=1:r
            candsCoords(j,1:2,i) = cands(j).Lmax;
            candsCoords(j,3:4,i) = cands(j).Bkg1;
            candsCoords(j,5:6,i) = cands(j).Bkg2;
            candsCoords(j,7:8,i) = cands(j).Bkg3;
        end
    elseif i>=100
        load([candsDirectory sprintf('\\cands%d.mat',i)])
        r = candsPerFrame(i);
        for j=1:r
            candsCoords(j,1:2,i) = cands(j).Lmax;
            candsCoords(j,3:4,i) = cands(j).Bkg1;
            candsCoords(j,5:6,i) = cands(j).Bkg2;
            candsCoords(j,7:8,i) = cands(j).Bkg3;
        end
    else
        disp('function cannot process >= 1000 files');
    end
end
if ~isempty(candsPerFrame)
    avg = mean(candsPerFrame);
    sigma = std(candsPerFrame);
end

