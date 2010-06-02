function [results] = stageFitLifetimes(data, histName, cutoff_start)
% STAGEFITLIFETIMES(DATA) fits the lifetimes contained in the structure in
% several stages, using data with different framerates. Globally averaged
% and leave-one-out histograms are returned for later jackknife estimation.
%
% SYNOPSIS [results] = stageFitLifetimes(data)
%
% INPUT     data:   structure containing all raw data
%
% OUTPUT:   results: structure containing the fields, for each framerate
%                       .hist_Xs
%                       .numcells_Xs
%                    where X is the framerate
%           lftHist_Xs :  (2+length(data)) x N matrices containing
%                                   row 1: time vector
%                                   row 2: normalized lifetimes
%                                   rows 3-end: leave-one-out normalized lifetimes
%
% Dinah Loerke, 31-Jul-2007
% Francois Aguet, Feb 2010

if nargin<2 || isempty(histName)
    histName = 'lftHist';
end

if nargin<3

    %========================================================================
    % 1. Detection artifact based on average of all data
    %=========================================================================
    
    % The detection artifact consists only of a small number of frames;
    % 100 frames are sufficient for its estimation.
    cutoff_art = 100;
    tvec = 1:cutoff_art;
    
    histMatrix_all = zeros(length(data), cutoff_art);
    idx = 1:cutoff_art;
    for i=1:length(data)
        if ~isempty(data(i).(histName))
            % normalize
            histMatrix_all(i,:) = data(i).(histName)(idx)/sum(data(i).(histName)(idx));
        end
    end
    
    % combined/averaged histograms, normalized
    cHistogramVector = nanmean(histMatrix_all,1);
    cHistogramVector = cHistogramVector/nansum(cHistogramVector);
    
    
    % plot results
    figure;
    plot(tvec, cHistogramVector, 'k.'); hold on;
    xlabel('frames');
    ylabel('average frequency');
    axis([0 cutoff_art 0 1.1*max(cHistogramVector(:))]);
    
    
    % fit averaged results with multiple exponentials (k = 1)
    prmVect = [0  0.3 0.3 1  0.3 20 1];
    estVect = [1  1 1 0  1 1 0];
    prmVect = fitNWeibull(tvec, cHistogramVector, prmVect, estVect, 'PDF');
    
    % sort parameter vector according to population mean
    N = 2;
    [~, order] = sort(prmVect(3:3:end));
    idx = reshape(2:3*N+1, [3 N]);
    idx = reshape(idx(:,order), [1 3*N]);
    prmVect(2:end) = prmVect(idx);
    [~, W] = nWeibull(tvec, prmVect, 'PDF');
    
    % plot first component, which is the detection artifact
    plot(tvec, W(1,:), 'r-');
    title('detection artifact');
    
    
    %=====================
    % Detection artifact
    %=======================
    % A detection artifact in form of a large number of short tracks (~1-4 frames)
    % is part of the liftime distributions. A cutoff to exclude this part of the
    % histogram is defined below. The cutoff is in principle independent of frame rate.
    
    % distance at which the artifact exponential decreases past 0.1% of its amplitude at 1 frame
    % prmVect(3)*log((1-exp(-1/prmVect(3)))/0.001) + 1
    cutoff_start = ceil(-prmVect(3)*log(0.001)+1) + 1; %%%%%%%%%%%% + 1 remove
    fprintf('Detection artifact: cutoff at %d frames.\n', cutoff_start);
end
results.detectionCutoff = cutoff_start;


% In this implementation, since we want to average over the data without
% any breaks due to averaging over different numbers of histograms, the
% histograms are all cut to the same length. Thus, we need to know what's
% the maximum common time span for all the movies.

% ========================================================================
% Average all data
% =========================================================================

timespan = [data.framerate].*[data.movieLength];
nTracks = zeros(1, length(data));
frameRates = unique([data.framerate]);
histStruct = struct('t', [], 'hist', []);

for r = 1:length(frameRates)
    idx = find([data.framerate]==frameRates(r) & arrayfun(@(x)~isempty(x.(histName)), data));
    
    if ~isempty(idx)
        nMovies = length(idx);
        if frameRates(r) >= 1
            frString = [num2str(frameRates(r)) 's'];
        else
            frString = [num2str(frameRates(r)*1000) 'ms'];
        end
        
        
        % end cutoff (long lifetimes close to the length of the movie have noisy
        % frequencies because of bad statistics) - in percent, e.g. 0.1 for 10%
        cutpoint = round(floor(min(timespan(idx))/frameRates(r))*0.9); % cut at 90% of total length
        %cutpoint = min([data(idx).movieLength]); % no cutoff
        
        % time vector, cut off at beginning and end
        tvec_curr = frameRates(r)*(cutoff_start+1:cutpoint-1);
        
        for p = 1:nMovies
            i = idx(p);

            % Number of trajectories in movie
            nTracks(i) = sum(data(i).(histName)(cutoff_start:end));
            
            % correct histogram for movie length
            currHistCorr = data(i).(histName) .* lftHist_correctionVector(data(i).movieLength);
            currHistCorr = currHistCorr(cutoff_start+1:cutpoint-1);
            
            if frameRates(r)==2
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
            else
                % pre-normalize first by remaining sum
                currHistNorm = currHistCorr/sum(currHistCorr);
            end
            
            % enter interpolated/renormalized histogram into summarizing matrix
            histStruct(p).t = tvec_curr;
            histStruct(p).hist = currHistNorm;
        end
        
        histMatrix = NaN(nMovies, length(tvec_curr));
        for p = 1:nMovies
            histMatrix(p,:) = histStruct(p).hist;
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
        tvecComp = tvec_curr(1)+frameRates(r)*(1:length(histVectorAve));
        tvecRes = tvecComp(finitepos);
        % define appropriate restricted lft vector
        histVectorAve = histVectorAve(finitepos);
        
        % leave-one-out averages for later jackknife
        fieldName = ['hist_' frString];
        if nMovies > 1
            histVectorAveN1 = zeros(nMovies, length(finitepos));
            for p = 1:nMovies
                histVectorAveN1(p,:) = nanmean(histMatrix([1:p-1 p+1:nMovies], finitepos), 1);
            end
            results.(fieldName) = [round(tvecRes*10)/10; histVectorAve; histVectorAveN1];
        else
            results.(fieldName) = [round(tvecRes*10)/10; histVectorAve];
        end
        
        results.(['numtracks_' frString]) = sum(nTracks(idx));
        results.histMatrix = histMatrix;
        
        
        % Plot cumulative lifetime distributions
        figure
        for p = 1:nMovies
            plot(histStruct(p).t, cumsum(histStruct(p).hist), 'b.-'); hold on;
        end
        plot(tvecRes, cumsum(histVectorAve),'r.-');
        axis([0 tvec_curr(end) 0 1.05]);
        xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 14);
        title(['Mean cumulative histogram ' frString]);
        
        % Plot lifetime distributions
        figure;
        %for p = 1:nMovies
        %    plot(histStruct(p).t, histStruct(p).hist, 'b.-'); hold on;
        %end        
        plot(tvecRes, histVectorAve, 'r.-');
        axis([0 tvecRes(end) 0 1.05*max(histVectorAve)]);
        xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 14);
        title(['Mean histogram ' frString]);
        
        if frameRates(r) == 2
            % Determine fraction of permanent objects
            prmVect = [0 0.3 4 2  0.3 8 2  0.3 90 1];
            estVect = [1 1 1 1 1 1 1 1 1 1];
            prmVect = fitNWeibull(tvecRes, 1-cumsum(histVectorAve), prmVect, estVect, 'PDF');
            
            figure;
            plot(tvecRes, 1-cumsum(histVectorAve),'k.-');
            hold on;
            np = (length(prmVect)-1)/3;
            colorV = ones(np,3);
            colorV(:,1) = (0:1/np:1-1/np);
            set(gca, 'ColorOrder', hsv2rgb(colorV), 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
            [w W] = nWeibull(tvecRes , prmVect, 'PDF');
            plot(tvecRes, W);
            plot(tvecRes, w, 'r');
            
            title('Immobile fraction offset fit');
            fprintf('Immobile fraction offset: %f\n', prmVect(1));
        end
    end
end