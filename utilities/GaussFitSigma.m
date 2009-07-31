function [sigma, sigmaList, fitImages] = GaussFitSigma(signalCell,sigmaGuess,fitOptions)
%GAUSSFITSIGMA fits the width of a 3D-Gaussian-like signal
%
% SYNOPSIS: [sigma, correction] = GaussFitSigma(signalCell,sigmaGuess,fitOptions)
%
% INPUT signalCell: cell array with signals. A signal is a 3D array that
%                   has the same size in x and y, and that contains exactly
%                   one Gaussian-like signal. The center of the array is
%                   assumed to be the center of the signal (interpolate
%                   from image if necessary - even with a rough initial
%                   guess for sigma you should be able to find the centroid
%                   quite reliably)
%		sigmaGuess: (opt) 1-by-2 array with initial guess for sxy and sz
%                   (in pixels). Default: [1,1]
%		fitOptions  (opt) Options for lsqnonlin. It is not recommended
%                   that you set 'Jacobian' to 'off'. See GaussFitND for
%                   details.
%
% OUTPUT sigma: 1-by-2 array with fitted widths of the Gaussian [sxy,sz].
%               parameters... outputs from GaussFitND. Needed to make
%               cdFitGauss easier to program.
%        sigmaList: list of individual sigmas
%        fitImages: nSignals-by-3 array with {raw signal, averaged signal, Gauss}
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 17-Sep-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================
%% CHECK INPUT
%=========================

% check signalCell
if nargin == 0 || isempty(signalCell) || iscell(signalCell) && any(cellfun('isempty',signalCell))
    error('GaussFitSigma needs nonempty signalCell')
end
% check for cell
if isnumeric(signalCell)
    tmp{1} = signalCell;
    signalCell = tmp;
    clear tmp
end
% count cells
nSignals = length(signalCell);

if nargout > 2
    collectFits = true;
    fitImages = cell(nSignals,3);
else
    collectFits = false;
end

% check for number of elements of signals. If all equal, we can set coords
% here. Otherwise, in loop
signalSize = cellfun('prodofsize',signalCell);
if all(signalSize(1)==signalSize)
    halfSize = (size(signalCell{1})-1)/2;
    [cx,cy,cz] = ndgrid(-halfSize(1):halfSize(1),...
        -halfSize(2):halfSize(2),-halfSize(3):halfSize(3));
    coordList = [cx(:),cy(:),cz(:)];
    recalcCoordList = false;
else
    recalcCoordList = true;
end

% check sigmaGuess
if nargin < 2 || isempty(sigmaGuess)
    sigmaGuess = [1 1];
else
    sigmaGuess = sigmaGuess(:)';
end

% check options
if nargin < 3
    fitOptions = [];
end

%=========================


%=========================
%% LOOP TO FIT SIGMA
%=========================

sigmaList = NaN(nSignals,2);
goodSignalCell = cell(size(signalCell));

for s=1:nSignals

    % read, average, norm
    % We don't care about noise-determined uncertainties here - all we want
    % is a good estimate for the width of the PSF. Since the sigma of the
    % Gauss is heavily correlated with the amplitude and the background,
    % these two parameters have to be fixed. The background is set to 0,
    % and the amplitude is chosen so that the integral of the Gaussian is
    % 1 (it is, of course, not exactly 1 in a finite volume, and the norm
    % is theoretically changing with sigma, though that error is
    % negligeable).
    % The background is determined by the wireframe of the signal array. In
    % order to avoid bias because of outliers or interfering signals, the
    % robustMean of the wireframe is taken as background, and outliers are
    % masked with NaN.
    % To increase SNR and to take care of uneven background, the signal is
    % averaged with rotations of itself. This assumes that the center of
    % the signal really is the center of the PSF, and that the PSF is
    % symmetric. If either assumption is violated, the sigma will get
    % larger. However, the effect is small for small violations, and in
    % general, overestimation of the size of the PSF leads to more stable
    % results (it may not allow to resolve a few overlapping tags in MMF,
    % but there won't be erroneous overfitting).


    % read signal
    signal = signalCell{s};

    % find wireframe - better: take outermost pixel layer. At 4 sigma,
    % there should be no intensity left above background, anyway.
    edgeMask = true(size(signal));
    edgeMask(2:end-1,2:end-1,2:end-1) = false;
    %     edgeMask(2:end-1,2:end-1,[1,end]) = false;
    %     edgeMask(2:end-1,[1,end],2:end-1) = false;
    %     edgeMask([1,end],2:end-1,2:end-1) = false;

    % estimate background
    edgeList = signal(edgeMask);
    % don't use robustMean. If there is intensity left over from a
    % centromere, it will be spread around with the averaging. Thus,
    % robustMean would underestimate the effective background
    %[background,dummy,dummy,outlierIdx] = robustMean(edgeList(isfinite(edgeList)));
    background = nanmean(edgeList);

    % remove outliers - useless, b/c outside the outermost layer, we won't
    % remove anything!
    %     edgeIdx = find(edgeList);
    %     signal(edgeIdx(outlierIdx)) = NaN;
    signal = signal-background;

    if collectFits
        % collect raw-bg
        fitImages{s,1} = signal./nansum(signal(:));
    end

    % average all (nanmean in case there were outliers etc.
    % average Z
    signal = nanmean(cat(4,signal,signal(:,:,end:-1:1)),4);
    % average Y
    signal = nanmean(cat(4,signal,signal(:,end:-1:1,:)),4);
    % average X
    signal = nanmean(cat(4,signal,signal(end:-1:1,:,:)),4);
    % average X/Y (rotation by 90°)
    signal = nanmean(cat(4,signal,permute(signal,[2,1,3])),4);

    % normalize intensity
    signal = signal./nansum(signal(:));

    if collectFits
        % collect averaged
        fitImages{s,2} = signal;
    end

    goodSignalCell{s} = signal;

end % loop preparing

% fit
for s=nSignals:-1:1
    % fit - request all output so that we can give output if necessary
    try
    if recalcCoordList
        halfSize = (size(signalCell{s})-1)/2;
        [cx,cy,cz] = ndgrid(-halfSize(1):halfSize(1),...
            -halfSize(2):halfSize(2),-halfSize(3):halfSize(3));
        coordList = [cx(:),cy(:),cz(:)];
    end
    [parameters,sigmaParameters,Q] = ...
        GaussFitND(goodSignalCell{s}(:), ...
        coordList, {'sxy', 's3'}, ...
        [0,0,0,NaN,sigmaGuess,0], 1, fitOptions);

    % relevant parameters are #5 and #6
    sigmaList(s,:) = parameters(1,5:6);
    catch
        sigmaList(s,:) = [];
        goodSignalCell(s) = [];
        signalCell(s) = [];
    end

end

%debug
if 0
    rows = ceil(sqrt(nSignals));
    cols = ceil(nSignals/rows);
    figure
    for s=1:nSignals
        subplot(rows,cols,s),
        imshow(cat(1,nanmean(signalCell{s},3),squeeze(nanmean(signalCell{s},1))'),[])
    end
    % xy, xz projections are blue/yellow - white means that there is
    % perfect overlap. Red-green is not so good for seeing the exact
    % overlap, but better for the faint differences at low signals.
    figure
    for s=1:nSignals
        subplot(rows,cols,s),
        gauss=GaussMask3D(sigmaList(s,[1,1,2]),size(goodSignalCell{s}),...
            [0,0,0],2);
        red = nanmean(goodSignalCell{s},3);
        maxRed = max(red(:));
        red = red/maxRed;
        green = red;
        blue = nanmean(gauss,3)/maxRed;
        imshow(cat(3,red,green,blue));
    end
    figure
    for s=1:nSignals
        subplot(rows,cols,s),
        gauss=GaussMask3D(sigmaList(s,[1,1,2]),size(goodSignalCell{s}),...
            [0,0,0],2);
        red = squeeze(nanmean(goodSignalCell{s},1))';
        maxRed = max(red(:));
        red = red/maxRed;
        green = red;
        blue = squeeze(nanmean(gauss,1))'/maxRed;
        imshow(cat(3,red,green,blue));
    end
end
%

% collect Gauss
if collectFits
    for s=1:nSignals
        fitImages{s,3}=GaussMask3D(sigmaList(s,[1,1,2]),size(goodSignalCell{s}),...
            [0,0,0],2);
    end
end



% use robustMean for sigma
if nSignals > 2
    sigma = robustMean(sigmaList,1);
elseif nSignals == 2
    sigma = mean(sigmaList,1);
else
    sigma = sigmaList;
end

% % debug
% fitSignal = mean(cat(4,goodSignalCell{:}),4);
% %
% % subplot(rows,cols,rows*cols)
% % imshow(mean(fitSignal,3),[])
% [parameters,sigmaParameters,Q] = ...
%     GaussFitND(fitSignal(:), coordList, {'sxy', 's3'}, [0,0,0,NaN,sigmaGuess,0], 1);
%
% % relevant parameters are #5 and #6
% sigmaFit = parameters(1,5:6);