function [mergedHistRes] = mergeFastSlowHistogramsPlat(Results, restrict, shapevec, estartvec, efixvec)
% merge the fast and slow histograms preserving the normalization, and
% fit multiple populations to lifetimes
% SYNOPSIS [mergedHistRes]=mergeFastSlowHistogramsPlat(Results, restrict, shape);
%
% INPUT     Results     = results of function StageFitLifetimePlat, which
%                         contain the summarized fast and slow results
%           restrict    = (optional) time restriction in seconds, e.g. 300
%           shape   = (optional) shape for populations, e.g. [2 2 1], where
%                   1 indicates an exponential distribution, and 
%                   2 indicates a Rayleigh distribution
%           estartvec (optional) = start values for fitting
%                   [ b a1 tau1 k1 a2 tau2 k2 .... an taun kn]
%           efixvec (optional) = vector with 1/0 values to indicate that
%                   certain start values will be fixed during the final
%                   fit; this vector needs to have the same length as
%                   startpar (3 times the number of populations plus one)
%                   and set to zero for no fixing, one for fixing, e.g.
%                   [ 1 0 0 0 0 0 0 ]
%                   for a fit with two populations where b is fixed, or
%                   [ 0 0 1 0 ]
%                   for a fit with one population where tau1 is fixed                  
%
% OUTPUT    mergedHistRes = merged Histogram results, which have the fields
%           .numcells   = number of trajectories, fast+slow;
%           .tvec_or    = time vector for original merged histogram
%           .hvec_or    =  original merged histogram
%           .tvec_cum   = time vector for original cumulative histogram
%           .hvec_cum   = cumulative merged histogram
%           .compactFitRes = compact matrix representation of fit results,
%               where the first column represents the relative
%               contributions of all populations (including the persistent
%               population offset), the second column represents the time 
%               constants, and the third column represents tau50   
%
%
% Dinah Loerke
% last modified August 20, 2008


% the input results have the following format:
tvecfast = Results.hist_400ms(1,:);
hvecfast = Results.hist_400ms(2,:);
tvecslow = Results.hist_2s(1,:);
hvecslow = Results.hist_2s(2,:);

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length

if nargin < 2 || isempty(restrict)
    restrict = max(tvecslow);
end

if nargin<3
    shapevec = [2 1 1];
end


% number of cells for statistics
nc_fast = Results.numtracks_400ms;
nc_slow = Results.numtracks_2s;

% offset from 1 in the cumulative histogram of the slow data
OffsetCH = 1 - sum(hvecslow);



% ========================================================================
% 
%    fit fast lifetime data to determine time constant of first Rayleigh
%
% =========================================================================

figure
plot(tvecfast, hvecfast,'mo');
axis([0.5 35 0 1.05*max(hvecfast)]); 

% fit histogram with 3 distributions, of which 2 are free, and the first
% one is kept fixed to a Rayleigh
startvF = [0    0.1 4   2   0.2 10  2   0.2 90  1];
fixvF   = [0    0   0   1   0   0   0   0   0   0]; 
estFast = fitcurveMultiWeibullODF_lsq( tvecfast, hvecfast, startvF, fixvF);

% set all except offest to positive
estFast(2:length(estFast)) = abs(estFast(2:length(estFast)));

for t = 1:3
    textF{t} = num2str(estFast(3*t-1:3*t+1));
end
text(10, 0.9*max(hvecfast), textF);
title('fast framerate');

% extract value of sigma for first fast Rayleigh distribution to be kept fixed in the slow data
sigRay = estFast(3);


% variation 03/21
% newly initialized startvector
%for s = 1:3
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
% area within the interval
norm_fast   = sum(hvecfast((tvecfast<=a2) & (tvecfast>a1)));
hnorm_fast  = hvecfast/norm_fast; % original histogram normalized with area
% framerate
dtfast = max(diff(tvecfast));

% SLOW
% positions within the specified time interval
% area within the interval
norm_slow  = sum(hvecslow((tvecslow<=a2) & (tvecslow>a1)));
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
tfrebin(tfcrop>=clo) = NaN;
hfrebin = hfcrop;
hfrebin(tfcrop>=clo) = NaN;
step = round(dtslow/dtfast);

pstart = find(tfcrop>=clo, 1, 'first');
ti = 0;
for t = pstart:step:(length(tfcrop)-1)
    tfrebin(pstart+ti) = tfcrop(t);
    hfrebin(pstart+ti) = nanmean(hfcrop(t:t+step-1));
    ti = ti+1;
end



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
%sfinal = [tvecslow(length(tvecslow)) hvecslow(length(hvecslow)) OffsetCH];
%snew = [tcomb(length(tcomb)) hcomb(length(hcomb)) OffsetNorm];


% renormalize the merged histogram including the newly normalized offset;
% the resulting cumulative histogram should again converge towards 1 or a
% slightly lower (but not higher!) value
% NOTE: for both the norm and for the cumulative distribution, take the
% different time steps in tcomb into account!!
dtvec = round(10*diff(tcomb))/10; 
dtvec = [dtvec(1) dtvec]; 

finalsum = sum(hcomb.*dtvec) + OffsetNorm*dtvec(length(dtvec));
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


% % plot the results
% figure
% plot(tfrebin,hfrebin/dtfast,'r.');
% hold on;
% plot(tscrop,hscrop/dtslow,'b.');
% plot(tcomb, hcomb, 'go');
% %axis([0 120 0 0.2]);
% 
% % figure;
% % subplot(2,1,1);
% % plot(tcomb,hfitNorm);
% % axis([0 maxt -0.001 0.03]);
% % subplot(2,1,2);
% % plot(tcomb,hfitNormCum);
% % axis([0 maxt -0.05 1.05]);
% figure;

% ========================================================================
% 
%       fit combined original histogram with with multiple distributions,
%       number and shape is determined by shape input vector
%
% =========================================================================


startv1template = [0  0.2 sigRay 2  0.1 10 2  0.2 90 1  0.2 100 1];
startv1 = startv1template(1:1+3*length(shapevec));
startv1(4:3:end) = shapevec;

fixv1template   = [1  0 1 1  0 0 1  0 0 1  0 0 1]; 
fixv1 = fixv1template(1:length(startv1));

[est1] = fitcurveMultiWeibullODF_lsq(tcomb, hfitNorm, startv1, fixv1);

% for t=1:length(shapevec), text1{t}=num2str(est1(3*t-1:3*t+1)); end
% text( 10, 0.9*max(hfitNorm),text1 );
% axis([0 maxt 0 1.01*max(hfitNorm)]);

% ========================================================================
% 
%       fit combined cumulative histogram with multiple distributions,
%       number and shape is determined by shape input vector
%
% =========================================================================


% figure
startv2 = est1;
startv2(4:3:length(startv2)) = shapevec;
% new addition: setting the fixvector position 1 to value 2 ensures that
% this parameter (the offset) is kept to a positive value during the
% fitting - this is owed to the fact that the data are normalized and the
% sum of the three individual populations cannot be more than 1 , only less
% than 1
fixv2template = [0  0  1  1  0  0  1  0  0  1  0  0  1]; 
fixv2 = fixv1;
fixv2(1:length(fixv2)) = fixv2template(1:length(fixv2));

% subplot(2,1,2); hold on;
axis([0 maxt 0 1.01*max(hfitNormCum)]);
est2 = fitcurveMultiWeibullCDF_lsq(tcomb, hfitNormCum, startv2, fixv2);

% for t=1:length(shapevec), text2{t}=num2str( round(1000*est2(3*t-1:3*t+1))/1000 ); end
% text( 10, 0.9*max(hfitNormCum),text2 );
% axis([0 maxt 0 1.01*max(hfitNormCum)]);
% title('merged cumulative histogram');

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

%figure
startv3 = est2; 
startv3(4:3:length(startv3)) = shapevec;
startv3(1) = 0;
fixv3template   = [0  0  0  1  0  0  1  0  0  1  0  0  1]; 
fixv3 = fixv2;
fixv3(1:length(fixv3)) = fixv3template(1:length(fixv3));

if nargin>3
    startv3 = estartvec;
    if nargin >4
        fixv3 = efixvec;
    end
end


axis([0 maxt 0 1.01*max(1-hfitNormCum)]);

[est3] = fitcurveMultiWeibullCDFi_lsq(tcomb, 1-hfitNormCum, startv3, fixv3);


% relative contributions
atpos = [1, 2:3:length(est3)];
ATcont = abs(est3(atpos))/sum(abs(est3(atpos)));
ATcont = round(ATcont*10000)/100;

% time constants
ATtau = est3(3:3:length(est3));
ATtau = [0 round(ATtau*100)/100];



for t = 1:length(shapevec)+1
    text3{t} = num2str([ATcont(t) ATtau(t)]);
end
text( 10, 0.9*max(1-hfitNormCum), text3 );
axis([0 maxt 0 1.01*max(1-hfitNormCum)]);


% box-and whisker ranges for the different distributions - the first
% distribution is the offset and therefore excluded
[range50(1)] = 0;
for p = 2:length(ATtau)
    [range50(p)] = boxwhiskerPerRange(ATtau(p), shapevec(p-1), 0.5);
    [range75(p)] = boxwhiskerPerRange(ATtau(p), shapevec(p-1), 0.75);
end
range50 = round(range50*100)/100;
range75 = round(range75*100)/100;


%OUTPUT
mergedHistRes.numcells = nc_fast + nc_slow;
mergedHistRes.tvec_or = tcomb;
mergedHistRes.hvec_or = hfitNorm;
mergedHistRes.tvec_cum = tcomb;
mergedHistRes.hvec_cum = hfitNormCum;

% Output display vector. Format: relative contributions p1 p2 p3 p4
format bank
mergedHistRes.compactFitRes = [ATcont' ATtau' range50' range75'];

%tplot = [0.5:0.5:max(t)];
%if min(nonzeros(diff(t(:))))<1
%    tplot = t;
%end

% ========================================================================
%                               display final results
% ========================================================================


figure;
% plot original histogram
plot(tcomb, hfitNorm, 'b.-');
axis([0 150 -0.001 1.05*max(hfitNorm(:))]);
title('merged orginal histogram');
xlabel('lifetime');
ylabel('frequency');
legendstring(1) = {'measured data'};

hold on;
% final fit results
estimates = est3;
p1params = [0 estimates(2:4)];
p1curve = multiWeibullODF(tcomb,p1params);
plot(tcomb,p1curve,'c-'); 
legendstring(2) = {'population 1'};
if length(estimates)>6
    p2params = [0 estimates(5:7)];
    p2curve = multiWeibullODF(tcomb,p2params);
    plot(tcomb,p2curve,'g-');
    legendstring(3) = {'population 2'};
    if length(estimates)>9
        p3params = [0 estimates(8:10)];
        p3curve = multiWeibullODF(tcomb,p3params);
        plot(tcomb,p3curve,'m-');
        legendstring(4) = {'population 3'};
        if length(estimates)>12
            p4params = [0 estimates(11:13)];
            p4curve = multiWeibullODF(tcomb,p4params);
            plot(tcomb,p4curve,'r-');
            legendstring(5) = {'population 4'};
        end
            
    end
end
legend(legendstring);