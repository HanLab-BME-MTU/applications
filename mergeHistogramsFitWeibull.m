function [mergedHistRes] = mergeHistogramsFitWeibull(Results, restrict, shape, estartvec, efixvec)
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
% Dinah Loerke
% last modified May 14, 2009


% loop over fast and slow results (multiple entries, which are the results
% of deleting one movie after the other)

tvecslow = Results.hist_2s(1,:);
tvecfast = Results.hist_400ms(1,:);

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length
if nargin<2
    restrict = max(tvecslow);
end
if nargin<3
    shape = [2 1 1];
end

% number of cells for statistics
nc_fast = Results.numtracks_400ms;
nc_slow = Results.numtracks_2s;

nFast = size(Results.hist_400ms,1)-1;
nSlow = size(Results.hist_2s,1)-1;

nc = 1;
mergedHistRes(1:(nFast+nSlow-1)) = struct();
% fast framerate: global average followed by leave-one-out averages
figure;
for nf = 1:nFast
    
    hvecfast = Results.hist_400ms(nf+1,:);
    
    if nf==1
        nslowloop = nSlow;
    else
        nslowloop = 1;
    end
    for ns = 1:nslowloop
        
        hvecslow = Results.hist_2s(ns+1,:);
        
        % offset from 1 in the cumulative histogram of the slow data
        OffsetCH = 1 - sum(hvecslow);
        
        
        % determine Rayleigh time constant from fast lifetime histogram, then
        % combine slow and fast and fit together
        
        
        % ========================================================================
        %       fit fast lifetime data to determine Rayleigh
        % =========================================================================
        
        plot(tvecfast, hvecfast, 'mo');
        axis([0.5 35 0 1.05*max(hvecfast)]);
        
        % fit histogram with 3 distributions, of which 2 are free, and the first
        % one is kept fixed to a Rayleigh
        startvF = [0  0.1 4 2  0.2 15 2  0.2 60 1];
        fixvF   = [0  0 0 1  0 0 0  0 0 0];
        %estFast = fitcurveMultiWeibullODF_lsq( tvecfast, hvecfast, startvF, fixvF);
        %sum((nWeibull(tvecfast, estFast, 'PDF')-hvecfast).^2)
        estFast = fitNWeibull(tvecfast, hvecfast, startvF, 1-fixvF, 'PDF');
        %sum((nWeibull(tvecfast, estFast, 'PDF')-hvecfast).^2)
        
        % figure
        % plot(tvecfast, hvecfast, 'k.-');
        % hold on;
        % plot(tvecfast, nWeibull(tvecfast, estFast1, 'PDF'), 'b');
        % plot(tvecfast, nWeibull(tvecfast, estFast2, 'PDF'), 'r');
        
        
        
        % set all except offest to positive
        estFast(2:end) = abs(estFast(2:end));
        
        text1 = cell(1,length(shape));
        for t = 1:length(shape)
            text1{t} = num2str(estFast(3*t-1:3*t+1), '%11.2f');
        end
        text(10, 0.9*max(hvecfast), text1);
        title('fast framerate');
        
        % extract value of sigma for first fast Rayleigh distribution to be kept fixed in the slow data
        sigRay = estFast(3);
        
        % ========================================================================
        %       make continuous distribution of fast and slow
        % =========================================================================
        
        % for this purpose, normalize both fast and slow histogram to the area
        % between a1 and a2 (in seconds)
        a1 = 28;
        a2 = 48;
        
        % FAST
        % area within the interval
        norm_fast   = sum(hvecfast((tvecfast<=a2) & (tvecfast>a1)));
        hnorm_fast  = hvecfast/norm_fast; % original histogram normalized with area
        % framerate
        dtfast = max(diff(tvecfast));
        
        % SLOW
        norm_slow  = sum(hvecslow((tvecslow<=a2) & (tvecslow>a1)));
        hnorm_slow  = hvecslow/norm_slow;
        % framerate
        dtslow = max(diff(tvecslow));
        
        
        % for data fitting, cut off last part of fast data and first part of slow data and combine them together
        
        % lower croplimit (in seconds)
        clo = a1;
        % higher crop limit (in seconds)
        chi = 100;
        
        hfcrop = hnorm_fast(tvecfast<=chi);
        tfcrop = tvecfast(tvecfast<=chi);
        
        hscrop = hnorm_slow(tvecslow>=clo);
        tscrop = tvecslow(tvecslow>=clo);
        
        
        % re-bin the portion of the fast histogram above clo to match the slow histogram
        tfrebin = tfcrop;
        tfrebin(tfcrop>=clo) = NaN;
        hfrebin = hfcrop;
        hfrebin(tfcrop>=clo) = NaN;
        step = round(dtslow/dtfast);
        
        pstart = find(tfcrop>=clo, 1, 'first');
        ti = 0;
        for t=pstart:step:(length(tfcrop)-1)
            tfrebin(pstart+ti) = tfcrop(t);
            hfrebin(pstart+ti) = nanmean(hfcrop(t:t+step-1));
            ti=ti+1;
        end
        
        plot(tfrebin,hfrebin/dtfast,'r.');
        hold on
        plot(tscrop,hscrop/dtslow,'b.');
        axis([0 120 0 0.2]);
        
        
        
        % === combined vectors   -  average
        f1 = find(tfrebin == clo);
        f2 = find(tfrebin == (chi-dtslow));
        s1 = find(tscrop == clo);
        s2 = find(tscrop == (chi-dtslow));
        % first stretch: taken from fast
        tcomb(1:f1) = tfrebin(1:f1);
        hcomb(1:f1) = hfrebin(1:f1)/dtfast;
        % second stretch: average of fast and slow
        tcomb(f1:f2) = tfrebin(f1:f2);
        % average is weighted with number of points in each category
        %wfast = nc_fast/(nc_fast+nc_slow);
        %wslow = nc_slow/(nc_fast+nc_slow);
        hcomb(f1:f2) = ( (nc_fast*(hfrebin(f1:f2)/dtfast)) + (nc_slow*(hscrop(s1:s2)/dtslow)))/(nc_fast+nc_slow);
        % third stretch: taken from slow
        par3 = hscrop(s2:length(hscrop))/dtslow; lp3 = length(par3);
        tcomb(f2:f2+lp3-1) = tscrop(s2:s2+lp3-1);
        hcomb(f2:f2+lp3-1) = hscrop(s2:s2+lp3-1)/dtslow;
        
        % normalized offset (taken from original slow histogram) to match the new
        % normalized values of the histogram: first divide by the same norm, and
        % then multiply by the framerate dtslow as done for the histogram above
        OffsetNorm = (OffsetCH/norm_slow)/dtslow;
        
        
        % renormalize the merged histogram including the newly normalized offset;
        % the resulting cumulative histogram should again converge towards 1 or a
        % slightly lower (but not higher!) value
        % NOTE: for both the norm and for the cumulative distribution, take the
        % different time steps in tcomb into account!!
        dtvec = round(10*diff(tcomb))/10;
        dtvec = [dtvec(1) dtvec];
        
        finalsum = sum(hcomb.*dtvec)+OffsetNorm*dtvec(end);
        hfitNorm = hcomb/finalsum;
        hfitNormCum = cumsum(hfitNorm.*dtvec);
        
        % restrict the vectors as specified
        pr = find(tcomb>restrict, 1, 'first');
        if ~isempty(pr)
            tcomb = tcomb(1:pr);
            hfitNorm = hfitNorm(1:pr);
            hfitNormCum = hfitNormCum(1:pr);
        end
        maxt = tcomb(end);
        
        plot(tcomb,hfitNorm);
        axis([0 maxt -0.001 0.03]);
        
        plot(tcomb,hfitNormCum);
        axis([0 maxt -0.05 1.05]);
        
        
        
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
        
        %est1 = fitcurveMultiWeibullODF_lsq(tcomb, hfitNorm, startv1, fixv1);
        %sum((hfitNorm - nWeibull(tcomb, est1, 'PDF')).^2)
        est1 = fitNWeibull(tcomb, hfitNorm, startv1, 1-fixv1, 'PDF');
        %sum((hfitNorm - nWeibull(tcomb, est1, 'PDF')).^2)
        
        % display
        for t = 1:length(shape)
            text1{t} = num2str(est1(3*t-1:3*t+1), '%11.2f');
        end
        text(10, 0.9*max(hfitNorm), text1);
        axis([0 maxt 0 1.01*max(hfitNorm)]);
        
        
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
        
        axis([0 maxt 0 1.01*max(hfitNormCum)]);
        %est2 = fitcurveMultiWeibullCDF_lsq( tcomb, hfitNormCum, startv2, fixv2);
        %sum((hfitNormCum - nWeibull(tcomb, est2, 'CDF')).^2)
        est2 = fitNWeibull(tcomb, hfitNormCum, startv2, 1-fixv2, 'CDF');
        %sum((hfitNormCum - nWeibull(tcomb, est2, 'CDF')).^2)
        
        for t = 1:length(shape)
            text1{t} = num2str(est2(3*t-1:3*t+1), '%11.2f');
        end
        text(10, 0.9*max(hfitNormCum), text1);
        axis([0 maxt 0 1.01*max(hfitNormCum)]);
        
        
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
        
        axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
        
        %est3 = fitcurveMultiWeibullCDFi_lsq(tcomb, 1-hfitNormCum, startv3, fixv3);
        %sum((1-hfitNormCum - nWeibull(tcomb, est3, 'SVF')).^2)
        est3 = fitNWeibull(tcomb, 1-hfitNormCum, startv3, 1-fixv3, 'SVF');
        %sum((1-hfitNormCum - nWeibull(tcomb, est3, 'SVF')).^2)
        
        
        % relative contributions
        ATcont = abs(est3([1 2:3:end]));
        ATcont = 100*ATcont/sum(ATcont);
        
        % time constants
        ATtau = [0 est3(3:3:end)];
        
        for t = 1:length(shape)+1
            text1{t} = num2str([ATcont(t) ATtau(t)], '%11.2f');
        end
        text(10, 0.9*max(1-hfitNormCum), text1);
        axis([0 maxt 0 1.01*max(1-hfitNormCum)]);
        
        % box-and whisker ranges for distributions (1st distribution is offset, excluded)
        range50 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.5), 2:length(ATtau))];
        range75 = [0 arrayfun(@(p) boxwhiskerPerRange(ATtau(p), shape(p-1), 0.75), 2:length(ATtau))];
        

        
        %OUTPUT
        mergedHistRes(nc).numcells = nc_fast + nc_slow;
        mergedHistRes(nc).tvec_or = tcomb;
        mergedHistRes(nc).hvec_or = hfitNorm;
        mergedHistRes(nc).tvec_cum = tcomb;
        mergedHistRes(nc).hvec_cum = hfitNormCum;
        mergedHistRes(nc).prmVect = est3;
        
        mergedHistRes(nc).p25 = ATtau(2:end) .* nthroot(-log(0.75), shape);
        mergedHistRes(nc).p75 = ATtau(2:end) .* nthroot(-log(0.25), shape);
        
        
        % Output display vector. Format: relative contributions p1 p2 p3 p4
        mergedHistRes(nc).compactFitRes = [ATcont' ATtau' range50' range75'];
        nc = nc+1;
    end
end

% ========================================================================
%                               display final results
% ========================================================================


colorV = [0 1 1; 0 1 0; 1 0 1; 1 0 0];
t = mergedHistRes(1).tvec_or;

[w W] = nWeibull(t, [0 mergedHistRes(1).prmVect(2:end)]);

sum((w-mergedHistRes(1).hvec_or).^2)

lg = arrayfun(@(x) ['Population ' num2str(x)], 1:size(W,1), 'UniformOutput', false);
lg = ['Measured data', 'Model', lg];

figure;
plot(mergedHistRes(1).tvec_or, mergedHistRes(1).hvec_or, 'b.-');
axis([0 150 -0.001 1.05*max(mergedHistRes(1).hvec_or)]);
set(gca, 'ColorOrder', colorV);

hold on;
plot(t, w, 'r');
plot(t, W);
legend(lg);
title('Merged histograms');
xlabel('Lifetime (s)');
ylabel('Relative frequency');
