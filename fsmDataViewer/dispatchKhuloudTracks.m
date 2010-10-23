function dataLayer = dispatchKhuloudTracks(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

load(fileList{1});

if ~exist('tracksFinal', 'var')
    error('Unable to find track feature info.');
end

trackSEL = getTrackSEL(tracksFinal);

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    % Find which tracks live in iFrame
    inFrameInd = find(trackSEL(:,1) <= iFrame & trackSEL(:,2) >= iFrame);
    
    tracks = cell(numel(inFrameInd),1);
    
    indX = 8*(iFrame-trackSEL(inFrameInd,1)) + 1;
    indY = 8*(iFrame-trackSEL(inFrameInd,1)) + 2;
    
    for iiTrack = 1:numel(inFrameInd)
        iTrack = inFrameInd(iiTrack);
       
        x = tracksFinal(iTrack).tracksCoordAmpCG(1:8:indX(iiTrack));
        y = tracksFinal(iTrack).tracksCoordAmpCG(2:8:indY(iiTrack));
        
        tracks{iiTrack} = [x', y'];
    end
    
    dataLayer{iFrame} = tracks;
end
