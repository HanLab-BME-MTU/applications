function idlisttrack = retrofitErrors(idlisttrack,idlist)
%RETROFITERRORS replaces tracker errors with detector errors
%
% SYNOPSIS [idlisttrack] = retrofitErrors(idlisttrack)
%
% INPUT    idlisttrack
%          
%
% OUTPUT   idlisttrack
%
% retrofitErrors replaces tracker errors with detector errors
% This is needed because we realized that the tracker errors were
% unreasonably small.
% For tags found by the detector we use their
% values found in detectQ_Pix and substitute them in totalQ_Pix.
% For tags not found by the detector we use an interpolated value
%
%
% Author: Eugenio 4/18/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTimePoints = length(idlisttrack);

% We first find the number of tags from a non-deleted frame
firstFrame = 0;
nFound = 0;
while ~nFound
    firstFrame = firstFrame +1;
    if ~isempty(idlist(1,firstFrame).linklist)
        nTags = length(idlisttrack(1,firstFrame).linklist(:,1));
        nFound = 1;
    end
end
% We first find the list of detected tags
listDetected = logical(zeros(3*nTags,nTimePoints));

for indexTag = 1:nTags
    for indexTimePoints = 1:nTimePoints
        % column 3 in linklist has detected spot labels 0:=detected,
        % 1:=missed
        % For convenience we expand it
        if ~isempty(idlist(1,indexTimePoints).linklist)
            listDetected(3*indexTag-2,indexTimePoints) = ~idlist(1,indexTimePoints).linklist(indexTag,3);
            listDetected(3*indexTag-1,indexTimePoints) = ~idlist(1,indexTimePoints).linklist(indexTag,3);
            listDetected(3*indexTag,indexTimePoints) = ~idlist(1,indexTimePoints).linklist(indexTag,3);
        end
    end
end

% We now get the mean Q-values for the detected tags
qValuesDetected = NaN(3*nTags,nTimePoints);

for indexTimePoints = 1:nTimePoints
    % column 3 in linklist has detected spot labels 0:=detected,
    % 1:=missed
    %         qValuesDetected(3*indexTag-2,indexTimePoints) = idlist(1,indexTimePoints).linklist(indexTag,3);
    %         qValuesDetected(3*indexTag-1,indexTimePoints) = ~idlist(1,indexTimePoints).linklist(indexTag,3);
    %         qValuesDetected(3*indexTag,indexTimePoints) = ~idlist(1,indexTimePoints).linklist(indexTag,3);
    % We check that we had a detected frame, with non-empty info and
    % linklist fields
    if ~isempty(idlist(1,indexTimePoints).info)&&~isempty(idlist(1,indexTimePoints).linklist)
        qValuesDetected(1:3*nTags,indexTimePoints) = diag(full(idlist(1,indexTimePoints).info.detectQ_Pix));
        if isfield(idlisttrack(1,indexTimePoints).info,'totalQ_Pix')
            qValuesTracked(1:3*nTags,indexTimePoints) = diag(full(idlisttrack(1,indexTimePoints).info.totalQ_Pix));
        end
    end
end

% We erase the values for tags not detected
qValuesDetected(~listDetected) = NaN;

% qValuesDetected increase gradually with time
% We interpolate the data with csaps
% plot(qValuesDetected')
splinedQValuesDetected = zeros(size(qValuesDetected));
warning('off','SPLINES:CHCKXYWP:NaNs')
for indexTag = 1:nTags
    splinedQValuesDetected(3*indexTag-2,1:nTimePoints) = csaps(1:nTimePoints,qValuesDetected(3*indexTag-2,:),0.8,1:nTimePoints);
    splinedQValuesDetected(3*indexTag-1,1:nTimePoints) = splinedQValuesDetected(3*indexTag-2,1:nTimePoints);
    splinedQValuesDetected(3*indexTag,1:nTimePoints) = csaps(1:nTimePoints,qValuesDetected(3*indexTag,:),0.8,1:nTimePoints);
end
warning('on','SPLINES:CHCKXYWP:NaNs')

% We replace the values for tags not detected with interpolated values
qValuesDetected(~listDetected) = splinedQValuesDetected(~listDetected);

for indexTimePoints = 1:nTimePoints
    if isfield(idlisttrack(1,indexTimePoints).info,'totalQ_Pix')
        idlisttrack(1,indexTimePoints).info.oldTotalQ_Pix = idlisttrack(1,indexTimePoints).info.totalQ_Pix;
        idlisttrack(1,indexTimePoints).info.totalQ_Pix = diag(qValuesDetected(:,indexTimePoints));
    end
end





aargh = 0;%











