function [results] = stageFitLifetimes(data, histName)
% STAGEFITLIFETIMES(DATA) fits the lifetimes contained in the structure in
% several stages, using data with slow and fast framerates. Globally averaged
% and leave-one-out histograms are returned for later jackknife estimation.
%
% SYNOPSIS [results] = stageFitLifetimes(data)
%
% INPUT     data:   structure containing all raw data
%
% OUTPUT:   results: structure containing the fields
%                       .hist_fast
%                       .hist_slow
%                       .numcells_slow
%                       .numcells_fast
%           lftHist_fast =  (2+length(data))xn matrix containing
%                           row 1: time vector
%                           row 2: normalized lifetimes
%                           rows 3-end: leave-one-out normalized lifetimes
%                           (for FAST acquisition)
%           lftHist_slow =  (2+length(data))xn matrix containing
%                           row 1: time vector
%                           row 2: normalized lifetimes
%                           rows 3-end: leave-one-out normalized lifetimes
%                           (for SLOW acquisition)
%
% last modified DATE: 31-Jul-2007 (Dinah)
% Francois Aguet, Feb 2010

if nargin < 2
    histName = 'lftHist';
end

%========================================================================
% 1. average all data to identify initial detection artifact
%=========================================================================

% NOTE: since we're interested in the initial decay of the detection
% artifact, it should be sufficient to fit only the 100 or so frames; the
% cutoff is defined as cutoff_art

cutoff_art = 100;
tvec = 1:cutoff_art;

histMatrix_all = zeros(length(data),cutoff_art);
idx = 1:cutoff_art;
for i=1:length(data)
    if ~isempty(data(i).(histName))
        % normalize over area betwen 0 and cutoff_art frames. FA: cmp normalization over full domain?
        histMatrix_all(i,:) = data(i).(histName)(idx)/sum(data(i).(histName)(idx));
    end
end

% combined/averaged histograms, normalized
cHistogramVector = nanmean(histMatrix_all,1);
cHistogramVector = cHistogramVector/nansum(cHistogramVector);


% plot results
figure;
plot(tvec, cHistogramVector, 'b.-'); hold on;
xlabel('frames');
ylabel('average frequency');
axis([0 cutoff_art 0 1.1*max(cHistogramVector(:))]);


% fit averaged results with multiple exponentials (weibulls with k==1)
guessvector = [0    0.3 0.3 1   0.3 5   1   0.3 20  1];
fixvector =   [0    0   0   1   0   0   1   0   0   1];
estAll = fitcurveMultiWeibullODF_lsq(tvec, cHistogramVector, guessvector, fixvector);
% take abs of all except offset
estAll(2:end) = abs(estAll(2:end));


% plot first component, which is the detection artifact
plot(tvec, abs(estAll(2))*(1/estAll(3))*exp(-tvec/estAll(3)), 'g-');
title('detection artifact');

% extract value of sigma for initial detection artifact
sigArt = estAll(3);
% distance at which the artifact exponential decreases past 0.1% of its amplitude at 1 frame
d01percent = -sigArt*log(0.001)+1;
disp(['0.1% distance: ', num2str(d01percent), 'frames']);


% ========================================================================
%                   define following data cutoffs
% ========================================================================

% initial cutoff (fast lifetimes that are cut off because they may
% represent noisy detection artifact) - in frames
cutoff_start = ceil(d01percent);
% NOTE: I allow the option of using DIFFERENT values for cutoff_start in
% the fast and slow movies further below

% end cutoff (long lifetimes close to the length of the movie have noisy
% frequencies because of bad statistics) - in percent, e.g. 0.1 for 10%
cutoff_slow = 0.1;

% ========================================================================
%       determine which entries in the data structure represent fast or
%       slow movies and store the positions
% =========================================================================

% In this implementation, since we want to average over the data without
% any breaks due to averaging over different numbers of histograms, the
% histograms are all cut to the same length. Thus, we need to know what's
% the maximum common time span for all the fast/slow movies - as we only
% input movies with the same framerate (2s slow, 0.4s fast), the knowledge
% of the last frame suffices

frThreshold = 1;
slowIdx = find([data.framerate] >  frThreshold & arrayfun(@(x)~isempty(x.(histName)), data));
fastIdx = find([data.framerate] <= frThreshold & arrayfun(@(x)~isempty(x.(histName)), data));
timespan = [data.framerate].*[data.movieLength];

mintimespan_slow = min(timespan(slowIdx));
mintimespan_fast = min(timespan(fastIdx));

% ========================================================================
% Average all fast/slow data
% =========================================================================

% the following procedure is performed twice, once for fast, once for slow
% movies; the respective interpolation time vectors were defined above

nTracks = zeros(1, length(data));

% =========================================================================
% Process fast data
% =========================================================================
if ~isempty(fastIdx)
    nMovies = length(fastIdx);
    ipFramerate = 0.4;
    maxl = floor(mintimespan_fast/0.4);
    
    cutpoint = round(maxl*0.9); % cut at 90% of total length
    histMatrix  = NaN(nMovies,1000);
    
    figure;
    for p = 1:nMovies
        i = fastIdx(p);

        % time vector, cut off at beginning and end
        tvec_curr = data(i).framerate*(cutoff_start+1:cutpoint-1);
                
        % Number of trajectories in movie
        nTracks(i) = sum(data(i).(histName)(cutoff_start:end));
        
        % correct histogram for movie length
        currCorrVec = lftHist_correctionVector(data(i).movieLength); % length: movielength-2
        %currHistCorr = data(i).(histName)(1:data(i).movieLength-2).*currCorrVec;
        currHistCorr = data(i).(histName)(cutoff_start+1:cutpoint-1).*currCorrVec(cutoff_start+1:cutpoint-1);

        % pre-normalize first by remaining sum
        currHistNorm = currHistCorr/sum(currHistCorr);
        
        % enter interpolated/renormalized histogram into summarizing matrix
        histMatrix(p,:) = NaN;
        histMatrix(p,1:length(currHistNorm)) = currHistNorm;
        
        % plot results
        plot(tvec_curr,cumsum(currHistNorm), 'b.-'); hold on;
        axis([0 120 0 1.05]);
    end
    
    % when the averaging takes place, that's a possible source for
    % artifacts for the first rayleigh population -
    % e.g. when the first defined point is averaged over 1 entry
    % and then the next one over 6 points - the first point is only if
    % the number of entries is >1 for more than 2 rows
    for ci = 1:20
        col = histMatrix(:,ci);
        n_nan = sum(isnan(col));
        n_def = nMovies - n_nan;
        if (n_nan>n_def) && (n_def>0)
            histMatrix(:,ci) = NaN;
        end
    end
    
    histVectorAve = nanmean(histMatrix,1);
    finitepos = find(isfinite(histVectorAve));
    
    % define appropriate restricted time vector
    tvecComp = tvec_curr(1)+ipFramerate*(1:length(histVectorAve));
    tvecRes = tvecComp(finitepos);
    % define appropriate restricted lft vector
    histVectorAve = histVectorAve(finitepos);
        
    % leave-one-out averages for later jackknife
    if nMovies > 1
        histVectorAveN1 = zeros(nMovies, length(finitepos));
        for p = 1:nMovies
            histVectorAveN1(p,:) = nanmean(histMatrix([1:p-1 p+1:nMovies], finitepos), 1);
        end
        results.hist_fast = [round(tvecRes*10)/10; histVectorAve; histVectorAveN1];
    else
        results.hist_fast = [round(tvecRes*10)/10; histVectorAve];
    end
    
    results.numcells_fast = sum(nTracks(fastIdx));
    plot(tvecRes, cumsum(histVectorAve),'r.-');
    title('Cumulative average histogram');
    
    figure;
    plot(tvecRes, histVectorAve,'r.-');
    axis([0 120 0 1.05*max(histVectorAve)]);
    title('Average histogram');
end

% =========================================================================
% Process slow data
% =========================================================================
if ~isempty(slowIdx)
    nMovies = length(slowIdx);
    ipFramerate = 2;
    maxl = floor(mintimespan_slow/2);
    
    cutpoint = round(maxl*(1-cutoff_slow)); % cut at 90%
    histMatrix  = NaN(nMovies,1000);
    
    figure;
    for p = 1:nMovies
        i = slowIdx(p);

        % time vector, cut off at beginning and end
        tvec_curr = data(i).framerate*(cutoff_start+1:cutpoint-1);

        % Number of trajectories in movie
        nTracks(i) = sum(data(i).(histName)(cutoff_start:end));

        % correct histogram for movie length
        currCorrVec = lftHist_correctionVector(data(i).movieLength); % length: movielength-2
        %currHistCorr = data(i).lftHist(1:data(i).movieLength-2).*currCorrVec;
        currHistCorr = data(i).(histName)(cutoff_start+1:cutpoint-1).*currCorrVec(cutoff_start+1:cutpoint-1);

        % number of persistent trajectories: trajectories longer than cutpoint
        if isfield(data(i), 'lftHistCut') && ~isempty(data(i).lftHistCut)
            numPers = sum(data(i).lftHistCut(cutpoint:end));
        else
            % threshold is in seconds, cutoff must be multiplied by framerate
            numPers = numPersistentField(data(i), data(i).framerate*cutpoint);
        end
        % normalize by the total sum, where the sum of the histogram plus
        % the number of persistent objects equals 100%
        nTracks(i) = nTracks(i) + numPers;

        % pre-normalize first by remaining sum
        currHistNorm = currHistCorr/(sum(currHistCorr)+numPers);
        
        % enter interpolated/renormalized histogram into summarizing matrix
        histMatrix(p,:) = NaN;
        histMatrix(p,1:length(currHistNorm)) = currHistNorm;
        
        % plot results
        plot(tvec_curr,cumsum(currHistNorm), 'b.-'); hold on;
        axis([0 600 0 1.05]);
    end
    
    % when the averaging takes place, that's a possible source for
    % artifacts for the first rayleigh population -
    % e.g. when the first defined point is averaged over 1 entry
    % and then the next one over 6 points - the first point is only if
    % the number of entries is >1 for more than 2 rows
    for ci = 1:20
        col = histMatrix(:,ci);
        n_nan = sum(isnan(col));
        n_def = nMovies - n_nan;
        if (n_nan>n_def) && (n_def>0)
            histMatrix(:,ci) = NaN;
        end
    end
    
    histVectorAve = nanmean(histMatrix,1);
    finitepos = find(isfinite(histVectorAve));
    
    % define appropriate restricted time vector
    tvecComp = tvec_curr(1)+ipFramerate*(1:length(histVectorAve));
    tvecRes = tvecComp(finitepos);
    % define appropriate restricted lft vector
    histVectorAve = histVectorAve(finitepos);
    
    % leave-one-out averages for later jackknife
    if nMovies > 1
        histVectorAveN1 = zeros(nMovies, length(finitepos));
        for p = 1:nMovies
            histVectorAveN1(p,:) = nanmean(histMatrix([1:p-1 p+1:nMovies], finitepos), 1);
        end
        results.hist_slow = [round(tvecRes*10)/10; histVectorAve; histVectorAveN1];
    else
        results.hist_slow = [round(tvecRes*10)/10; histVectorAve];
    end
    
    results.numcells_slow = sum(nTracks(slowIdx));
    plot(tvecRes, cumsum(histVectorAve),'r.-');
    title('Cumulative average histogram (slow)');
    
    figure;
    plot(tvecRes, 1-cumsum(histVectorAve),'r.-');
    guessvec = [0 0.3 4 2 0.3 8 2 0.3 90 1];
    fixvec = [0 0 0 0 0 0 0 0 0 0];
    estSlowInv = fitcurveMultiWeibullODF_lsq(tvecRes, 1-cumsum(histVectorAve), guessvec, fixvec);
    disp(['immobile fraction offset = ', num2str(estSlowInv(1))]);
    title('Population fit');
    
    figure;
    plot(tvecRes, histVectorAve,'r.-');
    axis([0 600 0 1.05*max(histVectorAve)]);
    title('Average histogram (slow)');
end