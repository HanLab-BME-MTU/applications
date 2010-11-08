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
% Francois Aguet, last modified November 5, 2010


% histogram field
hName = fieldnames(Results);
hName = hName{cellfun(@(x) ~isempty(x), regexp(hName, 'hist_'))};
suffix = hName(regexp(hName, 'hist', 'end')+1:end);

histM = Results.(hName);
tvec = histM(1,:);

%tvecslow = Results.hist_2s(1,:);
%tvecfast = Results.hist_400ms(1,:);

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
% nc_fast = Results.numtracks_400ms;
% nc_slow = Results.numtracks_2s;
nc = Results.(['numtracks' suffix]);


%nFast = size(Results.hist_400ms,1)-1;
%nSlow = size(Results.hist_2s,1)-1;
n = size(histM,1)-1; % number of histograms: global avg + l1o

nc = 1;
histRes(1:n-1) = struct();


% Global average followed by leave-one-out averages
figure;
for nf = 1:n
    
    hvec = histM(nf+1,:);
    % offset from 1 in the cumulative histogram of the slow data
    OffsetCH = 1 - sum(hvec);
    
    
    % determine Rayleigh time constant from fast lifetime histogram, then
    % combine slow and fast and fit together
    % ========================================================================
    %       fit fast lifetime data to determine Rayleigh
    % =========================================================================
    
    %plot(tvec, hvec, 'mo');
    %axis([0.5 35 0 1.05*max(hvec)]);
    
    % fit histogram with 3 distributions, of which 2 are free, and the first
    % one is kept fixed to a Rayleigh
    startvF = [0  0.1 4 2  0.2 15 2  0.2 60 1];
    fixvF   = [0  0 0 1  0 0 0  0 0 0];
    %estFast = fitcurveMultiWeibullODF_lsq( tvecfast, hvecfast, startvF, fixvF);
    %sum((nWeibull(tvecfast, estFast, 'PDF')-hvecfast).^2)
    estP1 = fitNWeibull(tvec, hvec, startvF, 1-fixvF, 'PDF');
    %sum((nWeibull(tvecfast, estFast, 'PDF')-hvecfast).^2)
    
%     figure
%     plot(tvec, hvec, 'k.-');
%     hold on;
%     [w W] = nWeibull(tvec, estP1, 'PDF');
%     plot(tvec, w, 'b');
%     hold on;
%     plot(tvec, W, 'b--');
%     title('P1 fit');
    
    % set all except offest to positive
    estP1(2:end) = abs(estP1(2:end));
    
%     text1 = cell(1,length(shape));
%     for t = 1:length(shape)
%         text1{t} = num2str(estP1(3*t-1:3*t+1), '%11.2f');
%     end
%     text(0.7*tvec(end), 0.9*max(hvec), text1);

    
    % extract value of sigma for first fast Rayleigh distribution to be kept fixed in the slow data
    sigRay = estP1(3);
    
    
    %         % ========================================================================
    %         %       make continuous distribution of fast and slow
    %         % =========================================================================
    %
    %         % for this purpose, normalize both fast and slow histogram to the area
    %         % between a1 and a2 (in seconds)
    %         a1 = 28;
    %         a2 = 48;
    %
    %         % FAST
    %         % area within the interval
    %         norm_fast   = sum(hvecfast((tvecfast<=a2) & (tvecfast>a1)));
    %         hnorm_fast  = hvecfast/norm_fast; % original histogram normalized with area
    %         % framerate
    %         dtfast = max(diff(tvecfast));
    %
    %         % SLOW
    %         norm_slow  = sum(hvecslow((tvecslow<=a2) & (tvecslow>a1)));
    %         hnorm_slow  = hvecslow/norm_slow;
    %         % framerate
    %         dtslow = max(diff(tvecslow));
    %
    %
    %         % for data fitting, cut off last part of fast data and first part of slow data and combine them together
    %
    %         % lower croplimit (in seconds)
    %         clo = a1;
    %         % higher crop limit (in seconds)
    %         chi = 100;
    %
    %         hfcrop = hnorm_fast(tvecfast<=chi);
    %         tfcrop = tvecfast(tvecfast<=chi);
    %
    %         hscrop = hnorm_slow(tvecslow>=clo);
    %         tscrop = tvecslow(tvecslow>=clo);
    %
    %
    %         % re-bin the portion of the fast histogram above clo to match the slow histogram
    %         tfrebin = tfcrop;
    %         tfrebin(tfcrop>=clo) = NaN;
    %         hfrebin = hfcrop;
    %         hfrebin(tfcrop>=clo) = NaN;
    %         step = round(dtslow/dtfast);
    %
    %         pstart = find(tfcrop>=clo, 1, 'first');
    %         ti = 0;
    %         for t=pstart:step:(length(tfcrop)-1)
    %             tfrebin(pstart+ti) = tfcrop(t);
    %             hfrebin(pstart+ti) = nanmean(hfcrop(t:t+step-1));
    %             ti=ti+1;
    %         end
    %
    %         plot(tfrebin,hfrebin/dtfast,'r.');
    %         hold on
    %         plot(tscrop,hscrop/dtslow,'b.');
    %         axis([0 120 0 0.2]);
    %
    %
    %
    %         % === combined vectors   -  average
    %         f1 = find(tfrebin == clo);
    %         f2 = find(tfrebin == (chi-dtslow));
    %         s1 = find(tscrop == clo);
    %         s2 = find(tscrop == (chi-dtslow));
    %         % first stretch: taken from fast
    %         tcomb(1:f1) = tfrebin(1:f1);
    %         hcomb(1:f1) = hfrebin(1:f1)/dtfast;
    %         % second stretch: average of fast and slow
    %         tcomb(f1:f2) = tfrebin(f1:f2);
    %         % average is weighted with number of points in each category
    %         %wfast = nc_fast/(nc_fast+nc_slow);
    %         %wslow = nc_slow/(nc_fast+nc_slow);
    %         hcomb(f1:f2) = ( (nc_fast*(hfrebin(f1:f2)/dtfast)) + (nc_slow*(hscrop(s1:s2)/dtslow)))/(nc_fast+nc_slow);
    %         % third stretch: taken from slow
    %         par3 = hscrop(s2:length(hscrop))/dtslow; lp3 = length(par3);
    %         tcomb(f2:f2+lp3-1) = tscrop(s2:s2+lp3-1);
    %         hcomb(f2:f2+lp3-1) = hscrop(s2:s2+lp3-1)/dtslow;
    %
    %         % normalized offset (taken from original slow histogram) to match the new
    %         % normalized values of the histogram: first divide by the same norm, and
    %         % then multiply by the framerate dtslow as done for the histogram above
    %         OffsetNorm = (OffsetCH/norm_slow)/dtslow;
    %
    %
    %         % renormalize the merged histogram including the newly normalized offset;
    %         % the resulting cumulative histogram should again converge towards 1 or a
    %         % slightly lower (but not higher!) value
    %         % NOTE: for both the norm and for the cumulative distribution, take the
    %         % different time steps in tcomb into account!!
    %         dtvec = round(10*diff(tcomb))/10;
    %         dtvec = [dtvec(1) dtvec];
    %
    %         finalsum = sum(hcomb.*dtvec)+OffsetNorm*dtvec(end);
    %         hfitNorm = hcomb/finalsum;
    %         hfitNormCum = cumsum(hfitNorm.*dtvec);
    %
    %         % restrict the vectors as specified
    %         pr = find(tcomb>restrict, 1, 'first');
    %         if ~isempty(pr)
    %             tcomb = tcomb(1:pr);
    %             hfitNorm = hfitNorm(1:pr);
    %             hfitNormCum = hfitNormCum(1:pr);
    %         end
    %         maxt = tcomb(end);
    %

    tcomb = tvec;
    
    %dtvec = diff(tvec);
    %dtvec = [dtvec(1) dtvec];
    
    dt = tvec(2)-tvec(1);
    %sum(hvec*dt)
    hfitNorm = hvec;
    hfitNormCum = cumsum(hfitNorm*dt);
    maxt = tvec(end);
    
    
    
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
    
    plot(tcomb, hfitNorm, 'k.-');
    hold on;
    plot(tvec, nWeibull(tvec, est1, 'PDF'), 'r');
    hold off;
%     title('est1 fit');
    
    
    % display
%     for t = 1:length(shape)
%         text1{t} = num2str(est1(3*t-1:3*t+1), '%11.2f');
%     end
%     text(0.7*tvec(end), 0.9*max(hfitNorm), text1);
    axis([0 maxt 0 1.05*max(hfitNorm)]);
    
    
%     % ========================================================================
%     %
%     %       fit combined cumulative histogram with multiple distributions,
%     %       number and shape is determined by shape input vector
%     %
%     % =========================================================================
%     
%     startv2 = est1;
%     startv2(4:3:end) = shape;
%     
%     fixv2template = [0  0 1 1  0 0 1  0 0 1  0 0 1  0 0 1];
%     fixv2 = fixv2template(1:length(fixv1));
%     
%     axis([0 maxt 0 1.01*max(hfitNormCum)]);
%     est2 = fitNWeibull(tcomb, hfitNormCum, startv2, 1-fixv2, 'CDF');
%     
%     figure
%     plot(tcomb, hfitNormCum, 'k.-');
%     hold on;
%     plot(tvec, nWeibull(tvec, est2, 'CDF'), 'r');
%     title('est2 fit');
%     
%     
%     for t = 1:length(shape)
%         text1{t} = num2str(est2(3*t-1:3*t+1), '%11.2f');
%     end
%     text(10, 0.9*max(hfitNormCum), text1);
%     axis([0 maxt 0 1.01*max(hfitNormCum)]);
    
    
%     % ========================================================================
%     %
%     %       fit inverse cumulative histogram with 3 distributions, with the
%     %       constraint that the fit can only have positive offset
%     %
%     % =========================================================================
%     
%     
%     % what constraints can we use?
%     % the constraint that the sum of the amplitudes (plus offset) in the
%     % cumulative histogram equals 100 is maybe not entirely justified, since
%     % the missing data points at the beginning of the fast acquisition can
%     % cause a shift. However, we can introduce the constraint that the limit
%     % (for infinity) of the sum of the distributions may be less, but not more
%     % than 100.
%     
%     startv3 = est2;
%     startv3(4:3:length(startv3)) = shape;
%     startv3(1) = 0;
%     fixv3template = [0  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1];
%     fixv3 = fixv3template(1:length(fixv2));
%     
%     if nargin>3
%         startv3 = estartvec;
%         if nargin>4
%             fixv3 = efixvec;
%         end
%     end
%     
%     axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
%     est3 = fitNWeibull(tcomb, 1-hfitNormCum, startv3, 1-fixv3, 'SVF');
%     
%     figure
%     plot(tcomb, hfitNorm, 'k.-');
%     hold on;
%     plot(tvec, nWeibull(tvec, est3, 'PDF'), 'r');
%     title('est3 fit');
    
    
    
    % relative contributions
    ATcont = abs(est1([1 2:3:end]));
    ATcont = 100*ATcont/sum(ATcont);
    
    % time constants
    ATtau = [0 est1(3:3:end)];
    
    %for t = 1:length(shape)+1
    %    text1{t} = num2str([ATcont(t) ATtau(t)], '%11.2f');
    %end
    %text(10, 0.9*max(1-hfitNormCum), text1);
    %axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
    
    % display
    for t = 1:length(shape)
        text1{t} = num2str(est1(3*t-1:3*t+1), '%11.2f');
    end
    text(0.7*tvec(end), 0.9*max(hfitNorm), text1);
    %axis([0 maxt 0 1.05*max(hfitNorm)]);
    
    
    
    % box-and whisker ranges for distributions (1st distribution is offset, excluded)
    range50 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.5), 2:length(ATtau))];
    range75 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.75), 2:length(ATtau))];
    
    
    
    %OUTPUT
    histRes(nc).numcells = nc;
    histRes(nc).tvec_or = tcomb;
    histRes(nc).hvec_or = hfitNorm;
    histRes(nc).tvec_cum = tcomb;
    histRes(nc).hvec_cum = hfitNormCum;
    histRes(nc).prmVect = est1;
    
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
