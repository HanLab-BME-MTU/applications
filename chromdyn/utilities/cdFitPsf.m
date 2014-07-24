function [fitStruct,fitImages] = cdFitPsf(fitStruct)
%CDFITPSF is a utility function that fits the width of a Gaussian to slist-spots
%
% in: fitStruct
%       .slist
%       .dataProperties
%       .movieDir
%       .rawMovieName
%       .spbIdx (optional; if omitted, all isolated spots will be fitted)
%       .waveIdx (optional, default = 1)
%
% out: fitStruct with updated sigmaCorrection and filter parameters in
%       fitStruct.dataProperties
%      fitImages: nSignals-by-3 list of raw, averaged and fitted images
%      (launches debug mode)

debug = false;
% check .spbIdx
if isfield(fitStruct,'spbIdx')
    checkSpb = true;
else
    checkSpb = false;
end
if ~isfield(fitStruct,'waveIdx') || isempty(fitStruct.waveIdx)
    fitStruct.waveIdx = 1;
end
if nargout > 1
    debug = true;
end

[dummy, dummy, loadStruct] = ...
    cdLoadMovie({fitStruct.rawMovieName,'corr/raw'}, ...
    fitStruct.movieDir, -1);
deltaFrames = loadStruct.frames2load{1}(1)-1;

% for isolated spots, position should always be ok (and, in
% principle, amplitude, though it will depend on sigma). To get
% correct sigma, fit width of isolated spots first
sigmaCorrectionOld = fitStruct.dataProperties.sigmaCorrection;

signalCell = {};

% loop through all frames. Check whether spbs are more
% than 7 sigmas away from everything else. If yes, load
% movie and read intensities from -4 sigma to +4 sigma
% No need to loop. A changed sigma may lead to a different
% number of spots being considered separated enough, but the
% fit seems to be stable and consistent already.
%             done = false;
originalSigma = fitStruct.dataProperties.FT_SIGMA([1,3])./sigmaCorrectionOld;
%             while ~done
for t = 1:length(fitStruct.slist)
    if ~isempty(fitStruct.slist(t).sp)
        slistCoords = cat(1,fitStruct.slist(t).sp.cord);
        
        % make distance matrix that counts in sigma. (metric
        % =1/sigma^2)
        metric = diag(1./fitStruct.dataProperties.FT_SIGMA(1:3).^2);
        spDist = distMat(slistCoords,metric);
        spDist = spDist + 10 * eye(size(spDist));
        
        % check closeness
        if checkSpb
            spList = cat(1,fitStruct.slist(t).sp.idxL);
            rowIdx = find(ismember(spList,fitStruct.spbIdx));
        else
            rowIdx = 1:length(spDist); %spDist is square
        end
        goodSpbIdx = find(all(spDist(rowIdx,:) > 7,2));
        
        if ~isempty(goodSpbIdx)
            
            % load movie
            cfLoadStruct = loadStruct;
            cfLoadStruct.frames2load{1} = t+deltaFrames;
            currentFrame = cdLoadMovie(cfLoadStruct.movieType, [], cfLoadStruct);
            currentFrame = currentFrame(:,:,:,fitStruct.waveIdx);
            
            % read intensities
            for s = 1:length(goodSpbIdx)
                spotCoord = slistCoords(rowIdx(goodSpbIdx(s)),[2,1,3]);
                hSize = ceil(fitStruct.dataProperties.FT_SIGMA*8)/2;
                [xx,yy,zz] = ndgrid(spotCoord(1)-hSize(1):spotCoord(1)+hSize(1),...
                    spotCoord(2)-hSize(2):spotCoord(2)+hSize(2),...
                    spotCoord(3)-hSize(3):spotCoord(3)+hSize(3));
                intList = interp3(currentFrame,yy,xx,zz,'*linear');
                
                % store in signalCell
                signalCell{end+1}=intList; %#ok<AGROW>
            end
        end
    end
end

% fit all
if ~isempty(signalCell)
    if debug
        [sigma,dummy,fitImages] = GaussFitSigma(signalCell,...
            fitStruct.dataProperties.FT_SIGMA([1,3]).*...
            fitStruct.dataProperties.sigmaCorrection);
    else
        sigma = GaussFitSigma(signalCell,...
            fitStruct.dataProperties.FT_SIGMA([1,3]).*...
            fitStruct.dataProperties.sigmaCorrection);
    end
    
    % update sigmaCorrection
    sigmaCorrection = ...
        sigma./originalSigma;
    
else
    sigmaCorrection = NaN(1,2);
    fitImages = [];
end


%disp(sigmaCorrection);
if isfinite(sigmaCorrection)
    fitStruct.dataProperties.sigmaCorrection = ...
        sigmaCorrection;
    fitStruct.dataProperties = defaultDataProperties(...
        fitStruct.dataProperties);
end