function [gaussFit,moreOutput] = cdFitGauss(inputStructure,fitList,allowOverlapForSigma)
%CDFITGAUSS fits Gaussians to idlist/slist data
%
% SYNOPSIS: [dataProperties,GaussWidth] = cdFitGauss(inputStructure)
%
% INPUT inputStructure: structure with fields
%		-idlist. Optional. If supplied, superfluous spots will be removed
%       -slist
%		-dataProperties
%       -movieDir
%		-rawMovieName
%       fitList: list of prarameters to fit (see GaussFitND)
%       allowOverlapForSigma : Normally, this is set to 0. Thus, the code
%           will not attempt to fit sigmas when the tags overlap.
%
% OUTPUT gaussFit: structure with the same size as slist with fields
%         - coords: fitParameters (x,y,z,a,sxy (or sx,sy),b)
%         - Q: cell array of covariance matrices (for every spot)
%         - chi2: ssq of residuals of fit
%        moreOutput: structure as gaussFit with field
%         - residualImage: image of residuals
%
% REMARKS gaussFit has the same spot order as slist, not idlist!!
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 06-Apr-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test inputStructure here

% test nargin
if nargin < 3 || isempty(allowOverlapForSigma)
    allowOverlapForSigma = false;
end

% update slist with spot order from idlist
if isfield(inputStructure,'idlist') && ~isempty(inputStructure.idlist)
    for t = 1:length(inputStructure.slist)
        if isempty(inputStructure.idlist(t).linklist)
            inputStructure.slist(t).sp = [];
        else
            % remove spots from slist that aren't interesting
            goodRows = inputStructure.idlist(t).linklist(:,2) > 0 & inputStructure.idlist(t).linklist(:,5) < 2;
            goodSpots = inputStructure.idlist(t).linklist(goodRows,2);
            rmIdx = setdiff(length(inputStructure.slist(t).sp),goodSpots);
            if ~isempty(rmIdx)
            inputStructure.slist(t).sp(rmIdx) = [];
            % remove spots from statistics
            % leave sigD for the moment
            inputStructure.slist(t).statistics.chi(rmIdx) = [];
            inputStructure.slist(t).statistics.snr(rmIdx) = [];
            inputStructure.slist(t).statistics.qAmp(rmIdx,:) = [];
            inputStructure.slist(t).statistics.qAmp(:,rmIdx) = [];
            idxList = (rmIdx(:)-1)*3;
            idxList = [idxList + 1; idxList + 2; idxList + 3];
            inputStructure.slist(t).statistics.Q(idxList,:) = [];
            inputStructure.slist(t).statistics.Q(:,idxList) = [];
            end
            
        end
    end
end

rawMovieName = fullfile(inputStructure.movieDir, inputStructure.rawMovieName);
dataProperties = inputStructure.dataProperties;
slist = inputStructure.slist;

% check fitList for 's' - if we fit sigmas, then we can't do overlapping
% spots
checkS = regexpi(fitList,'s');
if any([checkS{:}])
    sigmaFit = true;
else
    sigmaFit = false;
end
nFitParms = length(fitList);
checkB = regexpi(fitList,'b');
if any([checkB{:}])
    nFitParms = nFitParms - 1;
end
% nParameters is 3 coord, 1 amp, 2 sigma, 1 bg
nParameters = 7;


% decide on movieLoader
if isnumeric(inputStructure.rawMovieName)
    movieLoader = 'none';
elseif strcmpi(inputStructure.rawMovieName(end-2:end),'stk')
    movieLoader = 'Imaris';
elseif any(regexp(inputStructure.rawMovieName,'r3d$|3D.dv$'))
    movieLoader = 'cdLoadMovie';
else
    error('other movie loaders not implemented yet')
end



% load raw movie
switch movieLoader
    case 'cdLoadMovie'
        [rawMovie, movieHeader, loadStruct] = ...
            cdLoadMovie({rawMovieName,'corr/raw'}, [], dataProperties);
        % check if there are any leading darkframes we need to subtract
        deltaFrames = loadStruct.loadedFrames(1) - 1;
    case 'Imaris'
        if ~isfield(dataProperties,'crop')
            dataProperties.crop = [];
        end
        [rawMovie,movieSize,movieName,...
            moviePath,movieHeader,imarisHandle,loadStruct] = ...
            imarisImread(rawMovieName,[],dataProperties.crop,dataProperties.maxSize);
        deltaFrames = 0;
    case 'none'
        % read raw movie, set loadStruct and deltaFrames
        rawMovie = rawMovieName;
        loadStruct.loadedFrames = 1:size(rawMovie,5);
        loadStruct.frames2load = [];
        deltaFrames = 0;
        movieHeader.numTimepoints = size(rawMovie,5);
end

% loop through frames, then loop through spots in slist to fit sigmas

nTimepoints = length(slist);
gaussFit(1:nTimepoints) = struct('coords',[],'Q',{{}},'chi2',[]);
moreOutput(1:nTimepoints) = struct('residualImage',{{}});

for t=1:nTimepoints

    if ~isempty(slist(t).sp)
        % find current frame
        currentFrameIdx = find(t==loadStruct.loadedFrames-deltaFrames);
        if isempty(currentFrameIdx)
            % load next movie chunk
            switch movieLoader
                case 'cdLoadMovie'
                    [rawMovie, movieHeader, loadStruct] = ...
                        cdLoadMovie(loadStruct.movieType, [], loadStruct);
                case 'Imaris'
                    [rawMovie,dummy,dummy,...
                        dummy,dummy,dummy,loadStruct] = ...
                        imarisImread(loadStruct);
                case 'none'
                    error('more time points than movie frames!')
            end
            currentFrameIdx = find(t==loadStruct.loadedFrames-deltaFrames);
        end
        currentFrame = rawMovie(:,:,:,:,currentFrameIdx);

        % fit every spot individually. For a proper sigma-fit, spots should be
        % separated in space. Therefore, make sure that there isn't any spot too
        % close in slist

        slistCoords = [cat(1,slist(t).sp.cord),cat(1,slist(t).sp.amp),...
            cat(1,slist(t).sp.bg)];
        % !! flip x,y
        slistCoords = slistCoords(:,[2,1,3,4,5]);
        spotList = (1:size(slistCoords,1))';

        while ~isempty(slistCoords)
            % find overlapping spots
            [spotsIdx, mask] = discernspots(slistCoords(:,1:3),...
                dataProperties.movieSize(1:3),dataProperties);

            nSpots2fit = length(spotsIdx);

            if  ~allowOverlapForSigma && (nSpots2fit > 1 && sigmaFit)

                % don't fit - just write in NaN
                gaussFit(t).coords(spotList(spotsIdx),:) = NaN(nSpots2fit,nParameters);
                [gaussFit(t).Q{spotList(spotsIdx)}] = deal(NaN(nFitParms));
                if nargout > 1
                    [moreOutput(t).residualImage{spotList(spotsIdx)}]=...
                        deal(NaN);
                end
            else
                % fit multiple spots (or single spot)
                % read intensities, coordinates
                idxList=find(mask);
                mskData=currentFrame(idxList);
                [x,y,z] = ind2sub(dataProperties.movieSize(1:3),idxList);
                coordList = [x,y,z];

                % fit
                inputParameters = [slistCoords(spotsIdx,1:4),...
                    repmat(dataProperties.FT_SIGMA([1,3]),nSpots2fit,1),...
                    slistCoords(spotsIdx,5)];
                [parameters,sigmaParameters,Q,...
                    chiSquared,degreesOfFreedom,...
                    residualImage,resAndGauss] = ...
                    GaussFitND(mskData, coordList, ...
                    fitList, ...
                    inputParameters);

                % store output
                for i=1:nSpots2fit

                    % store output - mind the xy-switch!!
                    gaussFit(t).coords(spotList(spotsIdx(i)),:) = parameters(i,[2,1,3,4:nParameters]);
                    gaussFit(t).Q{spotList(spotsIdx(i))} = Q((i-1)*nFitParms+1:i*nFitParms,(i-1)*nFitParms+1:i*nFitParms);
                    if any(strcmpi(fitList,'x1')) && any(strcmpi(fitList,'x2'))
                    gaussFit(t).Q{spotList(spotsIdx(i))}(1:2,1:2) = ...
                        gaussFit(t).Q{spotList(spotsIdx(i))}([2,1],[2,1]);
                    end
                    gaussFit(t).chi2(spotList(spotsIdx(i)),1) = chiSquared;
                    % store residual image only if necessary
                    if nargout > 1
                        moreOutput(t).residualImage{spotList(spotsIdx(i))}=residualImage;
                    end

                end

                %             end
                %
                %         else
                %
                %             % read intensities, coordinates
                %             idxList=find(mask);
                %             mskData=currentFrame(idxList);
                %             [x,y,z] = ind2sub(dataProperties.movieSize(1:3),idxList);
                %             coordList = [x,y,z];
                %
                %             % fit
                %             inputParameters = [slistCoords(spotsIdx,1:4),...
                %                 dataProperties.FT_SIGMA([1,3]),...
                %                 slistCoords(spotsIdx,5)];
                %             [parameters,sigmaParameters,Q,...
                %                 chiSquared,degreesOfFreedom,...
                %                 residualImage,resAndGauss] = ...
                %                 GaussFitND(mskData, coordList, ...
                %                 fitList, ...
                %                 inputParameters);
                %
                %             disp(sprintf('%f %f %f %f %f %f %f',inputParameters))
                %             disp(sprintf('%f %f %f %f %f %f %f\n',parameters))
                %
                %             % store output
                %             gaussFit(t).coords(spotList(spotsIdx),:) = parameters;
                %             gaussFit(t).Q{spotList(spotsIdx)} = Q;
                %             gaussFit(t).chi2(spotList(spotsIdx),1) = chiSquared;
                %             % store residual image only if necessary
                %             if nargout > 1
                %                 moreOutput(t).residualImage{spotList(spotsIdx)}=residualImage;
                %             end

            end

            % remove fitted spots from list
            spotList(spotsIdx) = [];
            slistCoords(spotsIdx,:) = [];

        end % while there are coords to fit

        % write out new data
    end % if there is a non-empty slist

end % loop time

