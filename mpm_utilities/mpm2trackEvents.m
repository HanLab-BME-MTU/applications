function trackEvents = mpm2trackEvents(MPM,keepCompleteTracksOnly)

if nargin < 2 || isempty(keepCompleteTracksOnly)
    keepCompleteTracksOnly = true;
end

[nrows nFrames] = size(MPM);
nFrames = nFrames / 2;

%% Get the track mask
trackID = (1:nrows)';
trackMask = MPM(:,1:2:end) ~= 0;
    
if keepCompleteTracksOnly
    % remove any track that begins at 1st frame
    startAtFirstFrame = trackMask(:,1);
    trackMask(:,1) = false;
    for iFrame = 2:nFrames-1
        idxDead = trackID(~trackMask(:,iFrame));
        trackMask(:,iFrame) = trackMask(:,iFrame) & ~startAtFirstFrame;
        startAtFirstFrame(idxDead) = false;
    end
    % remove any track that ends at last frame
    endAtLastFrame = trackMask(:,end);
    trackMask(:,end) = false;
    for iFrame = nFrames-1:-1:1
        idxDead = trackID(~trackMask(:,iFrame));
        trackMask(:,iFrame) = trackMask(:,iFrame) & ~endAtLastFrame;
        endAtLastFrame(idxDead) = false;
    end
end

%% Get track events (i.e. firstFrame, lastFrame, row (in MPM))
trackEvents = cell(nFrames,1);

firstFrame = zeros(nrows,1);

for iFrame = 1:nFrames
    % keep first frame
    firstFrame(trackMask(:,iFrame) & firstFrame == 0) = iFrame;
 
    % Keep
    idxToKeep = ~trackMask(:,iFrame) & firstFrame ~= 0;
    
    trackEvents{iFrame} = [firstFrame(idxToKeep), repmat(iFrame-1,nnz(idxToKeep),1), trackID(idxToKeep)];
        
    % Reset
    firstFrame(~trackMask(:,iFrame)) = 0;
end

trackEvents = vertcat(trackEvents{:});
