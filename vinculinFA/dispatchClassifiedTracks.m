function dataLayer = dispatchClassifiedTracks(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

load(fileList{1});

if ~exist('tracksFinal', 'var')
    error('Unable to find track feature info.');
end

if ~exist('trackLabels', 'var')
    error('Unable to find track labels.');
end

nLabels = max(trackLabels);

colors = hsv(nLabels);

trackSEL = getTrackSEL(tracksFinal);

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    % Find which tracks live in iFrame
    inFrameInd = find(trackSEL(:,1) <= iFrame & trackSEL(:,2) >= iFrame);

    clear tracks;
    tracks(1:numel(inFrameInd)) = struct('tracksCoords', [], 'color',[]);
    
    indX = 8*(iFrame-trackSEL(inFrameInd,1)) + 1;
    indY = indX+1;
    
    for iiTrack = 1:numel(inFrameInd)
        iTrack = inFrameInd(iiTrack);
       
        x = tracksFinal(iTrack).tracksCoordAmpCG(1:8:indX(iiTrack));
        y = tracksFinal(iTrack).tracksCoordAmpCG(2:8:indY(iiTrack));

        isValid = ~isnan(x);
        
        tracks(iiTrack).trackCoords = [x(isValid)', y(isValid)'];
        tracks(iiTrack).color = colors(trackLabels(iTrack),:);
    end
    
    dataLayer{iFrame} = tracks;
end
