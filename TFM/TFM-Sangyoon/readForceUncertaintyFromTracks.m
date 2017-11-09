function tracksNA = readForceUncertaintyFromTracks(pathForMovieData,varargin)
% tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute) records
% pixel intensity from imgStack for entire track, even after ANA state in the last x,y
% position, and store it to attribute in tracksNA (1 means ampTotal, 2
% means forceMag
% assumes that tracks reflects coordinates from stage-drift-corrected
% images. imgStack is also from SDC output.
% get stack size
ip =inputParser;
ip.addParamValue('tracksNA',[],@isstruct); % selcted track ids
ip.parse(varargin{:});
tracksNA=ip.Results.tracksNA;

disp('Loading data ...')
tic
if isempty(tracksNA)
    tracksNA = load([pathForMovieData '/Colocalization/analysis1/data/tracksNA.mat'],'tracksNA');
    tracksNA = tracksNA.tracksNA;
end
tractionMaps = load([pathForMovieData '/TFMPackage/forceField/tractionMaps.mat'],'fCfdMap');
toc
fCfdMap = tractionMaps.fCfdMap;
curImg = fCfdMap;
numTracks = numel(tracksNA);
% parfor_progress(numel(tracksNA));
progressText(0,'Re-reading and tracking individual tracks:');

% parfor k=1:numel(tracksNA)
for k=1:numTracks
        try
            startFrame = tracksNA(k).startingFrameExtraExtra;
            endFrame = tracksNA(k).endingFrameExtraExtra;
        catch
            try
                startFrame = tracksNA(k).startingFrameExtra;
                endFrame = tracksNA(k).endingFrameExtra;
            catch
                startFrame = tracksNA(k).startingFrame;
                endFrame = tracksNA(k).endingFrame;
            end
        end
%         curRange = startFrame:endFrame;
%         tracksNA(k).forceUncertainty(curRange) = arrayfun(@(x) imgStack(round(tracksNA(k).yCoord(x)),round(tracksNA(k).xCoord(x)),x),curRange);
        for ii=startFrame:endFrame
            x = tracksNA(k).xCoord(ii);
            y = tracksNA(k).yCoord(ii);
            xi = round(x);
            yi = round(y);
            xRange = xi-1:xi+1;
            yRange = yi-1:yi+1;
            curAmpTotal = curImg(yRange,xRange);
            tracksNA(k).forceUncertainty(ii) = mean(curAmpTotal(:));
        end
    progressText(k/(numTracks-1),'Re-reading and tracking individual tracks:');
%     parfor_progress;
end
end
