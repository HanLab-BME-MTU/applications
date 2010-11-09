function [histRes] = fitWeibullToHist(Results, restrict, shape, estartvec, efixvec)
% merge the fast and slow histograms preserving the normalization, and
% fit multiple populations to lifetimes - PLUS estimation of experiment-
% to-experiment (i.e. cell-to-cell) error
% SYNOPSIS [mergedHistRes]=mergeFastSlowHistogramsPlat_E2E(Results, restrict, shape);
%
% INPUT     Results     = results of function StageFitLifetimePlat, which
%                         contain the summarized fast and slow results
%           restrict    = (optional) time restriction in seconds, e.g. 300
%           shape       = (optional) shape for populations, e.g. [2 2 1], where
%                       1 indicates an exponential distribution, and
%                       2 indicates a Rayleigh distribution
%           estartvec (optional) = start values for fitting
%                       [ b a1 tau1 k1 a2 tau2 k2 .... an taun kn]
%           efixvec (optional) = vector with 1/0 values to indicate that
%                       certain start values will be fixed during the final
%                       fit; this vector needs to have the same length as
%                       startpar (3 times the number of populations plus one)
%                       and set to zero for no fixing, one for fixing, e.g.
%                       [ 1 0 0 0 0 0 0 ]
%                       for a fit with two populations where b is fixed, or
%                       [ 0 0 1 0 ]
%                       for a fit with one population where tau1 is fixed
%
% OUTPUT    mergedHistRes   = merged Histogram results, which have the fields
%           .numcells       = number of trajectories, fast+slow;
%           .tvec_or        = time vector for original merged histogram
%           .hvec_or        =  original merged histogram
%           .tvec_cum       = time vector for original cumulative histogram
%           .hvec_cum       = cumulative merged histogram
%           .compactFitRes  = compact matrix representation of fit results,
%               where the first column represents the relative
%               contributions of all populations (including the persistent
%               population offset), the second column represents the time
%               constants, and the third column represents tau50
%
%
% Francois Aguet, last modified November 8, 2010


% histogram field
hName = fieldnames(Results);
hName = hName{cellfun(@(x) ~isempty(x), regexp(hName, 'hist_'))};
suffix = hName(regexp(hName, 'hist', 'end')+1:end);

histM = Results.(hName);
tvec = histM(1,:);

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length
if nargin<2
    restrict = max(tvec);
end
if nargin<3
    shape = [2 1 1];
end

% number of cells for statistics
nCells = Results.(['numtracks' suffix]);
nHist = size(histM,1)-1; % number of histograms: global avg + l1o

nc = 1;
histRes(1:nHist-1) = struct();

% histograms: global average followed by leave-one-out averages
figure;
for nf = 1:nHist
    
    hvec = histM(nf+1,:);
    
    % offset from 1 in the cumulative histogram
    OffsetCH = 1 - sum(hvec);
    
    % ========================================================================
    % Fit lifetime histogram to determine mean of 1st Rayleigh distribution
    % =========================================================================
    hold off;
    plot(tvec, hvec, 'k.-');
    hold on;
    %axis([0 35 0 1.05*max(hvec)]);
    
    % Fit histogram with 3 distributions (2 free, 1 fixed as Rayleigh)
    startvF = [0  0.1 4 2  0.2 15 2  0.2 60 1];
    fixvF   = [0  0 0 1  0 0 0  0 0 0];
    estFast = fitNWeibull(tvec, hvec, startvF, 1-fixvF, 'PDF');
       
    %text1 = cell(1,length(shape));
    %for t = 1:length(shape)
    %    text1{t} = num2str(estFast(3*t-1:3*t+1), '%11.2f');
    %end
    %text(0.7*tvec(end), 0.9*max(hvec), text1);
    
    % sigma of 1st Rayleigh distribution to be kept fixed
    sigRay = estFast(3);
    OffsetNorm = OffsetCH;
    
    % renormalize the merged histogram including the newly normalized offset;
    % the resulting cumulative histogram should again converge towards 1 or a
    % slightly lower (but not higher!) value
    % NOTE: for both the norm and for the cumulative distribution, take the
    % different time steps in tcomb into account!!
    dt = tvec(2)-tvec(1);
    
    hcomb = hvec;
    tcomb = tvec;
    
    finalsum = sum(hcomb*dt)+OffsetNorm*dt;
    hfitNorm = hcomb/finalsum;
    hfitNormCum = cumsum(hfitNorm*dt);
    
    % restrict the vectors as specified
    pr = find(tcomb>restrict, 1, 'first');
    if ~isempty(pr)
        tcomb = tcomb(1:pr);
        hfitNorm = hfitNorm(1:pr);
        hfitNormCum = hfitNormCum(1:pr);
    end
    maxt = tcomb(end);
    
    %plot(tcomb,hfitNorm);
    %axis([0 maxt -0.001 0.03]);
    
    %plot(tcomb,hfitNormCum);
    %axis([0 maxt -0.05 1.05]);
    
    
    % ========================================================================
    %
    %       fit combined original histogram with with multiple distributions,
    %       number and shape is determined by shape input vector
    %
    % =========================================================================
    
    
    
    startv1template = [0  0.2 sigRay 2  0.2 15 2  0.2 60 1  0.2 120 1  0.2 40 1];
    startv1 = startv1template(1:(1+3*length(shape)));
    startv1(4:3:end) = shape;
    
    fixv1template   = [1  0 1 1  0 0 1  0 0 1  0 0 1  0 0 1];
    fixv1 = fixv1template(1:length(startv1));    
    est1 = fitNWeibull(tcomb, hfitNorm, startv1, 1-fixv1, 'PDF');
    
    % display
    %for t = 1:length(shape)
    %    text1{t} = num2str(est1(3*t-1:3*t+1), '%11.2f');
    %end
    %text(0.7*tvec(end), 0.9*max(hfitNorm), text1);
    axis([0 tcomb(end) 0 1.01*max(hfitNorm)]);
    
    
    % ========================================================================
    %
    %       fit combined cumulative histogram with multiple distributions,
    %       number and shape is determined by shape input vector
    %
    % =========================================================================
    
    startv2 = est1;
    startv2(4:3:end) = shape;
    
    fixv2template = [0  0 1 1  0 0 1  0 0 1  0 0 1  0 0 1];
    fixv2 = fixv2template(1:length(fixv1));
    
    %axis([0 maxt 0 1.01*max(hfitNormCum)]);
    est2 = fitNWeibull(tcomb, hfitNormCum, startv2, 1-fixv2, 'CDF');
    
    %for t = 1:length(shape)
    %    text1{t} = num2str(est2(3*t-1:3*t+1), '%11.2f');
    %end
    %text(0.7*tvec(end), 0.9*max(hfitNormCum), text1);
    %axis([0 maxt 0 1.01*max(hfitNormCum)]);
    
    
    % ========================================================================
    %
    %       fit inverse cumulative histogram with 3 distributions, with the
    %       constraint that the fit can only have positive offset
    %
    % =========================================================================
    
    
    % what constraints can we use?
    % the constraint that the sum of the amplitudes (plus offset) in the
    % cumulative histogram equals 100 is maybe not entirely justified, since
    % the missing data points at the beginning of the fast acquisition can
    % cause a shift. However, we can introduce the constraint that the limit
    % (for infinity) of the sum of the distributions may be less, but not more
    % than 100.
    
    startv3 = est2;
    startv3(4:3:length(startv3)) = shape;
    startv3(1) = 0;
    fixv3template = [0  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1];
    fixv3 = fixv3template(1:length(fixv2));
    
    if nargin>3
        startv3 = estartvec;
        if nargin>4
            fixv3 = efixvec;
        end
    end
    
    %axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
    est3 = fitNWeibull(tcomb, 1-hfitNormCum, startv3, 1-fixv3, 'SVF');
        
    % relative contributions
    ATcont = abs(est3([1 2:3:end]));
    ATcont = 100*ATcont/sum(ATcont);
    
    % time constants
    ATtau = [0 est3(3:3:end)];
    
    for t = 1:length(shape)+1
        text1{t} = num2str([ATcont(t) ATtau(t)], '%11.2f');
    end
    text(0.7*tvec(end), 0.9*max(hfitNorm), text1);
    %axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
    
    % box-and whisker ranges for distributions (1st distribution is offset, excluded)
    range50 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.5), 2:length(ATtau))];
    range75 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.75), 2:length(ATtau))];
    
    
    
    %OUTPUT
    histRes(nc).numcells = nCells;
    histRes(nc).tvec_or = tcomb;
    histRes(nc).hvec_or = hfitNorm;
    histRes(nc).tvec_cum = tcomb;
    histRes(nc).hvec_cum = hfitNormCum;
    histRes(nc).prmVect = est3;
    
    histRes(nc).p25 = ATtau(2:end) .* nthroot(-log(0.75), shape);
    histRes(nc).p75 = ATtau(2:end) .* nthroot(-log(0.25), shape);
    
    
    % Output display vector. Format: relative contributions p1 p2 p3 p4
    histRes(nc).compactFitRes = [ATcont' ATtau' range50' range75'];
    nc = nc+1;
end


% ========================================================================
%                               display final results
% ========================================================================


colorV = [0 1 1; 0 1 0; 1 0 1; 1 0 0];
t = histRes(1).tvec_or;

[w W] = nWeibull(t, [0 histRes(1).prmVect(2:end)]);

sum((w-histRes(1).hvec_or).^2)

lg = arrayfun(@(x) ['Population ' num2str(x)], 1:size(W,1), 'UniformOutput', false);
lg = ['Measured data', 'Model', lg];

figure;
plot(histRes(1).tvec_or, histRes(1).hvec_or, 'b.-');
axis([0 150 -0.001 1.05*max(histRes(1).hvec_or)]);
set(gca, 'ColorOrder', colorV);

hold on;
plot(t, w, 'r');
plot(t, W);
legend(lg);
title('Merged histograms');
xlabel('Lifetime (s)');
ylabel('Relative frequency');
