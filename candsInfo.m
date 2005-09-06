function [candsPerFrame,avg,sigma] = candsInfo
%candsInfo queries user for cands directory and returns various stats
%
%SYNOPSIS   [candsPerFrame,avg,sigma] = candsInfo
%
%INPUT      none       
%
%OUTPUT     candsPerFrame: column vector, nth entry of which is the number
%                          of candidate speckles in file n
%           avg: mean number of candidate speckles per file, computed over 
%                n files
%           sigma: standard deviation of speckle number over the n frames
%
%c: 9-05 kathryn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[listOfFiles,directory] = searchFiles2('cands','.mat',[],'ask',0,'all');
sprintf('Directory searched:');
directory
sprintf('Number of cands files:');
numFiles = length(listOfFiles)

avg = 0;
avg = 0;
candsPerFrame = zeros(numFiles,1);


% put number of speckles per cands file in vector candsPerFrame
% ifelse used to accommodate cands file numbering system
for i=1:numFiles
    if i<10
        load([directory sprintf('\\cands00%d.mat',i)])
        candsPerFrame(i) = length(cands);
    elseif i>=10 & i<100
        load([directory sprintf('\\cands00%d.mat',i)])
        candsPerFrame(i) = length(cands);
    elseif i>=100
        load([directory sprintf('\\cands00%d.mat',i)])
        candsPerFrame(i) = length(cands);
    else
        disp('function cannot process >= 1000 files');
    end
end

avg = mean(candsPerFrame);
sigma = std(candsPerFrame);