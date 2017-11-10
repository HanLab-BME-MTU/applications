function [ tracksFinal ] = trackCloseGapsKalmanSparseSplitLAP(splitSize,movieInfo,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose)
%trackCloseGapsKalmanSparseSplitLAP takes runs trackCloseGapsKalmanSparse on
%pieces of movieInfo and then merges the resulting tracks 
%
% Input:
%       splitSize is the number frames per fraction       
%       
%       see trackCloseGapsKalmanSparse for a description of the above inputs
% except for splitSize
%
% Output:
%        tracksFinal, merges tracks
%
% Written by Jeffrey Werbin
% Harvard Medical School
% 11/10/13
%

% sets the division points for the movie
movieLen = numel(movieInfo);

if mod(movieLen,splitSize) > 0
    lim =[0:splitSize:movieLen,movieLen];
else
    lim = 0:splitSize:movieLen;
end

track = cell(numel(lim)-1,1);

for i = 1:numel(track)
    tic
    [temp,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(lim(i)+1:lim(i+1)),...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
    timeElap = toc;
    display(['For i = ',num2str(i),' it took ',num2str(timeElap),' seconds.']);
    track{i}=temp;    
end

tic
tracksFinal = track{1};
for i = 2:numel(track)
   tracksFinal = tracksFinalMergerLAP(tracksFinal,track{i},lim(i),costMatrices,gapCloseParam);    
end
timeElap = toc;
display(['Merging took ',num2str(timeElap),' seconds.']);

end

