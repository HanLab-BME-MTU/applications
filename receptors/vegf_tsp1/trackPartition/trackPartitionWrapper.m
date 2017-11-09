function trackPartitionWrapper(obj,varargin)
% Pass variables and parameters from TrackPartitionProcess to its
% functions, and output the partitioned tracks, diffusion analysis result
% (if the runDiffAnalsys parameter was set to true), and info on mean
% inside track length and number of inside tracks
%   [tracksPart,diffAnalysisRes,info] = trackPartitionWrapper(MD,params,trackChannel,maskChannel)

% Get input from Process
MD = obj.owner_;
maskMDPath = obj.maskMovieDataPath_;
maskMDFileName = obj.maskMovieDataFileName_;
params = obj.funParams_;
trackChannel = obj.trackChannel_;
maskChannel = obj.maskChannel_;
outFilePath = obj.outFilePaths_;

% Load mask MD
maskMD = load([maskMDPath,filesep,maskMDFileName]);
maskMD = maskMD.MD;

% Get MovieInfo and tracks
movieInfoFiles = maskMD.processes_{maskMD.getProcessIndex('SubResolutionProcess')}.outFilePaths_;
tracksFiles = MD.processes_{MD.getProcessIndex('TrackingProcess')}.outFilePaths_;
movieInfo = load(movieInfoFiles{maskChannel},'movieInfo');
movieInfo = movieInfo.movieInfo;
tracks = load(tracksFiles{trackChannel},'tracksFinal');
tracks = tracks.tracksFinal;

% If using existing mask
if params.useExistingMask
    % Check if one exists
    maskFile = MD.processes_{MD.getProcessIndex('TrackPartitionProcess')}.mask_;
    if ~isempty(maskFile)
        mask = load(maskFile);
        mask = mask.mask;
        % Check the upscale factor
        upscale = size(mask,1)/maskMD.imSize_(1);
        assert(mod(upscale,1) == 0,'The mask is not the same size as or a multiple of the size of the image')
        % Just upscale the tracks
        tracks = trackScalar(tracks,upscale);
    else
        % If none exists, create a mask and upscale tracks as normal
        [mask,tracks] = trackPartitionInit(tracks,movieInfo,maskMD,...
                 params.gaussianThresh,params.minMaskDiam,params.upscale);
    end
else
    % Create mask and upscale tracks
    [mask,tracks] = trackPartitionInit(tracks,movieInfo,maskMD,...
        params.gaussianThresh,params.minMaskDiam,params.upscale);
end

% Save mask
[maskPath,~,~] = fileparts(outFilePath{trackChannel});
maskPath = [maskPath,filesep,'partitionMask.mat'];
save(maskPath,'mask');
obj.mask_ = maskPath;

tracksPartTemp = trackPartition(tracks,mask,params.minTrackLength,...
    params.analysisStart,params.analysisEnd);

% Un-upscale tracks
tracksPart = trackScalar(tracksPartTemp,1/(params.upscale));

% Get the subset of inside tracks
subset = arrayfun(@(x) x.isInside,tracksPart);
info.tracksInside = tracksPart(subset);
info.nInsideTracks = sum(subset);

longest = 1;
totalLength = 0;
totalInsideTracks = 0;
for iTrack = 1:sum(subset)
    length = size(tracksPart(iTrack).tracksFeatIndxCG,2);
    if length > longest
        longest = length;
    end
    totalLength = totalLength+length;
    totalInsideTracks = totalInsideTracks+1;
end
meanLength = totalLength/totalInsideTracks;
info.meanLength = meanLength;

% Run diffusion analysis if needed
if params.runDiffAnalysis == true
    diffAnalysisRes = trackDiffusionAnalysis1(tracksPart); 
end

% Save
if params.runDiffAnalysis == true
    save(outFilePath{trackChannel},'tracksPart','info','diffAnalysisRes');
else
    save(outFilePath{trackChannel},'tracksPart','info');
end
end