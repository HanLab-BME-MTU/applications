function [mergedHistRes]=mergeFastSlowHistogramsPlat_E2E(Results, restrict, shape, estartvec, efixvec)
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

nfast = size(Results.hist_400ms, 1);
nslow = size(Results.hist_2s, 1);

nc = 1;

for nf = 1:nfast-1
    
    tvecfast = Results.hist_400ms(1,:);
    hvecfast = Results.hist_400ms(nf+1,:);
    
    if nf==1
        nslowloop = nslow-1;
    else
        nslowloop = 1;
    end
    
    for ns = 1: nslowloop
        tvecslow = Results.hist_2s(1,:);
        hvecslow = Results.hist_2s(ns+1,:);


        % optional input restrict: restrict final fitting analysis to a stretch of
        % the data, e.g. to the first 300s of the histogram, since measured
        % frequncies for lifetimes >300s are more than 50% speculative due to the
        % correction for movie length
        if nargin>1
            resT = restrict;
        else
            resT = max(tvecslow);
        end


        % number of cells for statistics
        nc_fast = Results.numtracks_2s;
        nc_slow = Results.numtracks_400ms;

        % offset from 1 in the cumulative histogram of the slow data
        OffsetCH = 1 - sum(hvecslow);


% determine Rayleigh time constant from fast lifetime histogram, then
% combine slow and fast and fit together


% ========================================================================
% 
%       fit fast lifetime data to determine Rayleigh
%
% =========================================================================


        % figure
        plot(tvecfast, hvecfast,'mo');
        axis([0.5 35 0 1.05*max(hvecfast)]); 

        % fit histogram with 3 distributions, of which 2 are free, and the first
        % one is kept fixed to a Rayleigh
        startvF = [0    0.1 4   2   0.2 15  2   0.2 60  1];
        fixvF   = [0    0   0   1   0   0   0   0   0   0]; 
        estFast = fitcurveMultiWeibullODF_lsq( tvecfast, hvecfast, startvF, fixvF);

        % set all except offest to positive
        estFast(2:length(estFast)) = abs(estFast(2:length(estFast)));

        for t=1:3, textF{t}=num2str(estFast(3*t-1:3*t+1)); end
        text( 10, 0.9*max(hvecfast), textF);

        title('fast framerate');

        % extract value of sigma for first fast Rayleigh distribution to be kept
        % fixed in the slow data
        sigRay = estFast(3);
        % sigRay = estFast(6);

        % variation 03/21
        % newly initialized startvector
        %for s=1:3
        %    shapePars(s) = min( max(round(estFast(1+s*3)),1) , 2);
        %end
        %niStartvector = estFast;
        %niStartvector(4:3:10) = shapePars;



% ========================================================================
% 
%       make continuous distribution of fast and slow
%
% =========================================================================

        % for this purpose, normalize both fast and slow histogram to the area
        % between a1 and a2 (in seconds)
        a1 = 28;
        a2 = 48;

        % FAST
        % positions within the specified time interval
        pnorm_fast  =  (tvecfast<=a2) & (tvecfast>a1) ;
        % area within the interval
        norm_fast   = sum(hvecfast(pnorm_fast));
        % original histogram normalized with area
        hnorm_fast  = hvecfast/norm_fast;
        % framerate
        dtfast = max(diff(tvecfast));

        % SLOW
        % positions within the specified time interval
        pnorm_slow  =  (tvecslow<=a2) & (tvecslow>a1) ;
        % area within the interval
        norm_slow  = sum(hvecslow(pnorm_slow));
        % original histogram normalized with area
        hnorm_slow  = hvecslow/norm_slow;
        % framerate
        dtslow = max(diff(tvecslow));


        % for data fitting, cut off last part of fast data and first part of slow
        % data and combine them together

        % lower croplimit (in seconds)
        clo = a1;
        % higher crop limit (in seconds)
        chi = 100;

        hfcrop = hnorm_fast(tvecfast<=chi);
        tfcrop = tvecfast(tvecfast<=chi);

        hscrop = hnorm_slow(tvecslow>=clo);
        tscrop = tvecslow(tvecslow>=clo);


        % re-bin the portion of the fast histogram above clo to match the slow
        % histogram
        tfrebin = tfcrop;
        tfrebin(tfcrop>=clo) = nan;
        hfrebin = hfcrop;
        hfrebin(tfcrop>=clo) = nan;
        step = round(dtslow/dtfast);

        pstart = find(tfcrop>=clo, 1, 'first');
        ti = 0;
        for t=pstart:step:(length(tfcrop)-1)
            tfrebin(pstart+ti) = tfcrop(t);
            hfrebin(pstart+ti) = nanmean(hfcrop(t:t+step-1));
            ti=ti+1;
        end


        % plot the results

        % figure
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

        hfit = hcomb;
        tfit = tcomb;


        % normalized offset (taken from original slow histogram) to match the new
        % normalized values of the histogram: first divide by the same norm, and 
        % then multiply by the framerate dtslow as done for the histogram above
        OffsetNorm = (OffsetCH/norm_slow)/dtslow;
        %sfinal = [tvecslow(length(tvecslow)) hvecslow(length(hvecslow)) OffsetCH];
        %snew = [tfit(length(tfit)) hfit(length(hfit)) OffsetNorm];


        % renormalize the merged histogram including the newly normalized offset;
        % the resulting cumulative histogram should again converge towards 1 or a
        % slightly lower (but not higher!) value
        % NOTE: for both the norm and for the cumulative distribution, take the
        % different time steps in tfit into account!!
        dtvec = round(10*diff(tfit))/10; 
        dtvec = [dtvec(1) dtvec]; 

        finalsum = sum(hfit.*dtvec)+OffsetNorm*dtvec(length(dtvec));
        hfitNorm = hfit/finalsum;
        hfitNormCum = cumsum(hfitNorm.*dtvec);




        % restrict the vectors as specified
        pr = find(tfit>resT, 1, 'first');
        if ~isempty(pr)
            tfit = tfit(1:pr);
            hfitNorm = hfitNorm(1:pr);
            hfitNormCum = hfitNormCum(1:pr);
        end



        maxt = tfit(length(tfit));
        % figure;
        %subplot(2,1,1);
        plot(tfit,hfitNorm);
        axis([0 maxt -0.001 0.03]);
        %subplot(2,1,2);
        plot(tfit,hfitNormCum);
        axis([0 maxt -0.05 1.05]);



        % ========================================================================
        % 
        %       fit combined original histogram with with multiple distributions,
        %       number and shape is determined by shape input vector
        %
        % =========================================================================

        if nargin>2
            shapevec = shape;
        else
            shapevec = [2 1 1];
        end

        startv1template = [0 0.2 sigRay 2 0.2 15 2 0.2 60 1 0.2 120 1 0.2 40 1];
        startv1 = zeros(1,1+3*length(shapevec));
        startv1(1:length(startv1)) = startv1template(1:length(startv1));
        startv1(4:3:length(startv1)) = shapevec;

        fixv1template   = [1 0 1 1 0 0 1 0 0 1 0 0 1 0 0 1]; 
        fixv1 = startv1;
        fixv1(1:length(fixv1)) = fixv1template(1:length(fixv1));

        %subplot(2,1,1); hold on;
        [est1] = fitcurveMultiWeibullODF_lsq( tfit, hfitNorm, startv1, fixv1);

        for t=1:length(shapevec), text1{t}=num2str(est1(3*t-1:3*t+1)); end
        text( 10, 0.9*max(hfitNorm),text1 );
        axis([0 maxt 0 1.01*max(hfitNorm)]);


        % ========================================================================
        % 
        %       fit combined cumulative histogram with multiple distributions,
        %       number and shape is determined by shape input vector
        %
        % =========================================================================


        hcfit = hfitNormCum;
        tcfit = tfit;


        % figure
        startv2 = est1;
        startv2(4:3:length(startv2)) = shapevec;
        % new addition: setting the fixvector position 1 to value 2 ensures that
        % this parameter (the offset) is kept to a positive value during the
        % fitting - this is owed to the fact that the data are normalized and the
        % sum of the three individual populations cannot be more than 1 , only less
        % than 1
        fixv2template   = [0 0 1 1 0 0 1 0 0 1 0 0 1 0 0 1]; 
        fixv2 = fixv1;
        fixv2(1:length(fixv2)) = fixv2template(1:length(fixv2));

        %subplot(2,1,2); hold on;
        axis([0 maxt 0 1.01*max(hcfit)]);
        [est2] = fitcurveMultiWeibullCDF_lsq( tcfit, hcfit, startv2, fixv2);

        for t=1:length(shapevec), text2{t}=num2str( round(1000*est2(3*t-1:3*t+1))/1000 ); end
        text( 10, 0.9*max(hcfit),text2 );
        axis([0 maxt 0 1.01*max(hcfit)]);


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

        % figure
        startv3 = est2; 
        startv3(4:3:length(startv3)) = shapevec;
        startv3(1) = 0;
        fixv3template   = [0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1]; 
        fixv3 = fixv2;
        fixv3(1:length(fixv3)) = fixv3template(1:length(fixv3));
        
        
        if nargin>3
            startv3 = estartvec;
            if nargin>4
                fixv3 = efixvec;
            end
        end

        
        axis([0 maxt 0 1.01*max(1-hcfit)]);
        
        [est3] = fitcurveMultiWeibullCDFi_lsq( tcfit, 1-hcfit, startv3, fixv3);


        % relative contributions
        atpos = [1, 2:3:length(est3)];
        ATcont = abs(est3(atpos))/sum(abs(est3(atpos)));
        ATcont = round(ATcont*1000)/10;

        % time constants
        ATtau = est3(3:3:length(est3));
        ATtau = [ 0 round(ATtau*10)/10];



        for t=1:length(shapevec)+1, text3{t}=num2str([ATcont(t) ATtau(t)]); end
        text( 10, 0.9*max(1-hcfit),text3 );
        axis([0 maxt 0 1.01*max(1-hcfit)]);


        % box-and whisker ranges for the different distributions - the first
        % distribution is the offset and therefore excluded
        [range50(1)] = 0;
        for p=2:length(ATtau)
            [range50(p)] = boxwhiskerPerRange(ATtau(p), shapevec(p-1), 0.5);
            [range75(p)] = boxwhiskerPerRange(ATtau(p), shapevec(p-1), 0.75);
        end
        range50 = round(range50*10)/10;
        range75 = round(range75*10)/10;


        %OUTPUT
        mergedHistRes(nc).numcells = nc_fast + nc_slow;
        mergedHistRes(nc).tvec_or = tfit;
        mergedHistRes(nc).hvec_or = hfitNorm;
        mergedHistRes(nc).tvec_cum = tcfit;
        mergedHistRes(nc).hvec_cum = hcfit;

        % output display vector
        % format: relative contributions p1 p2 p3 p4
        % compResMat = [ATcont' ATtau'];
        compResMat = [ATcont' ATtau' range50' range75'];
        format bank

        mergedHistRes(nc).compactFitRes = compResMat;
        nc = nc+1;
    
    end % of for ns
end % of for nf