function [BICres] = simpleBICtest(data)
% simpleBICtest performs a simple BIC test on the cumulative merged 
% histogram of the data specified
% SYNOPSIS [BICres] = simpleBICtest(data)
%
%INPUT      data:       movie data in the format prodcued by e.g.
%                       loadConditionData
%
%OUTPUT:    BICres:     results structure containing
%                       .startparams    = guess vector for fit;
%                       .fixparams      = fix vector 
%                       .estimates      = result of fit;
%                       .residuals      = residuals of fit
%                       .BICvalue       = Bayesian information value


% fill data structure to the required degree
data = determineMovieLength(data);
data = fillStructLifetimeHist(data);

% average and merge
res     = stageFitLifetimesPlat(data);
mer     = mergeFastSlowHistogramsPlat(res, 300, [2 2 1]);

% perform BIC test
% set start values for parameters
startpar = [0.1 0.2 10 2 0.2 20 2 0.2 90 1 0.2 25 1 0.2 35 1];

[BICres] = multiWeibullBICtest_cum(mer.tvec_cum,mer.hvec_cum,startpar);

end


%% ========================= SUBFUNCTION  =================================




function [res]=multiWeibullBICtest_cum(tvec,histvec,startvec)
%multiWeibullBICtest determines how many distributions are contained in the
%lifetime histogram histvec -  due to the peculiarities of my system, the
%first time constant is known and kept fixed
%SYNOPSIS [res]=multiWeibullBICtest(tvec,histvec,startvec)
%
%INPUT      tvec:       time vector
%           histvec:    corresponding lifetime vector (cumulative!)
%           startvec:   startvector for fitting
%
%OUTPUT:    res:        results structure containing
%                       .startparams    = guess vector for fit;
%                       .fixparams      = fix vector 
%                       .estimates      = result of fit;
%                       .residuals      = residuals of fit
%                       .BICvalue       = Bayesian information value


%% =======================================================================
%  for the final slow results, go through multiple steps to determine how
%  many distributions are needed for an optimum fit (based on Bayesian
%  information criterion)
% =======================================================================


% number of levels, i.e. distributions to test
numlev = 5;



%% kmat contains the shape parameter

% kmat contains as rows all the possible combinations of different
% distributions (k values) - identical permutations are already removed.
% The first distribution is kept fixed to the Rayleigh distribution 
% determined from the fast data, and subsequent distributions are added 
% that can be either Rayleighian (k=2) or exponential (k=1); or none (k=0)

kmat  = zeros(1,numlev);
kmat(1,1) = 1;
kmat(2,1) = 2;

% loop over number of distributions
for c=2:numlev
    [sx,sy] = size(kmat);
    % loop over (increasing) number of rows
    for r=1:sx
        %previous entry in this row
        prev_id = kmat(r,c-1);
        % if the previous distribution in this line is >0, add the
        % corresponding number of entries to the end of kmat
        if prev_id>0
           for a=1:prev_id
               [ssx,ssy] = size(kmat);
               kmat(ssx+1,:) = kmat(r,:);
               kmat(ssx+1,c) = a;
           end % of for a
        end % of if
    end % of for r
end % of for c


%% amat contains the values for the amplitude, lambdamat contains the
%% values for the lambda/sigma

% 1. amplitudes are split evenly between number of distributions
% 2. lambda for the first Rayleigh is known and kept fixed, for the other
% distributions, the additional populations are set in increments of 10
% frames


[sx,sy] = size(kmat);
for r=1:sx
    cnum = length(nonzeros(kmat(r,:)));
    amat(r,:) = 100*(1/cnum)*(kmat(r,:)>0);
    %if cnum>1
        lambdamat(r,1:cnum) = [10:20:(20*cnum-10)];
    %end % of if
end
%lambdamat(:,1) = sigRay;


%% for each row in the parameter matrices - which represents one specific
%% combination of distributions - perform a multiple-Weibull fit

BICtemp = 0;

h1 = waitbar(0,'fitting multiple distributions');
figure

for r=1:sx
    
    waitbar(r/sx);
    
    % to get the startvector and fixvector for this row, only use the
    % larger-than-zero values
    pos = find(kmat(r,:)>0);
    kvec = kmat(r,pos);
    avec = amat(r,pos);
    lambdavec = lambdamat(r,pos);
    startmat_curr = [avec; lambdavec; kvec];
    startvec_curr = [0 startmat_curr(:)'];
    if nargin>2
        ml = min(length(startvec),length(startvec_curr));
        startvec_curr(1:ml) = startvec(1:ml);
        startvec_curr(4:3:length(startvec_curr)) = kvec;
    end
    if length(startvec_curr)>=7
        if startvec_curr(7) == 2
            startvec_curr(6) = 2*startvec_curr(6);
        end
    end
    fixmat_curr = [0*avec; 0*lambdavec; (kvec>0)];
    fixvec_curr = [0 fixmat_curr(:)'];
    % NOTE: the sigma of the first Rayleigh (pos 3) is also kept fixed
    %fixvec_curr(3) = 1;
    
    axis([1 300 0 100]);
    % fit with multiple Weibulls; turn off display

    [estimates, residuals, estimatesSigma] = fitcurveMultiWeibullCDF_lsq(tvec, histvec, startvec_curr, fixvec_curr);
    %[estimates, residuals, estimatesSigma] = fitcurveMultiWeibullODF_lsq(tvec, histvec, startvec_curr, fixvec_curr);
    
    % NOTE: to prevent physcially meaningless populations, only positive
    % amplitudes and positive time constants are allowed
    estimates(2:length(estimates)) = abs(estimates(2:length(estimates)));
    
    fitresults(r).tvec = tvec;
    fitresults(r).hist = histvec;
    fitresults(r).startparams = startvec_curr;
    fitresults(r).fixparams = fixvec_curr;
    fitresults(r).estimates = estimates;
    fitresults(r).estimatesSig = estimatesSigma;
    fitresults(r).residuals = residuals;
    
    % value according to the Bayesian information criterion for this
    % combination
    
    % n = sample size
    BIC_n = length(find(isfinite(histvec)>0));
    % k = number of free parameters
    BIC_k = length(find(fixvec_curr==0));
    % rss = residual sum of squares
    BIC_rss = sum(residuals.^2);
    
    BICvalue = (BIC_n*log(BIC_rss/BIC_n)) + BIC_k*log(BIC_n);
    fitresults(r).BICvalue = BICvalue;
    
    % ====================================================================
    % =====     display results: comment/uncomment as desired
    
    % display values in separate figure only when BIC is a new minimum   
    if BICvalue<BICtemp
        
        
        plot(tvec, histvec,'bo');
        [wplot]=multiWeibullCDF(tvec,estimates);
        hold on; plot(tvec, wplot,'r-');
        titlestring=( [num2str(length(pos)),' populations, BIC=',num2str(BICvalue)] );
        title(titlestring);
        % relative amplitudes for text
        tamps = estimates([2:3:length(estimates)]);
        tramps= round(100*(tamps/sum(tamps)))/100;
        tsigs = round(100*estimates([3:3:length(estimates)]))/100;
        tshap = estimates([4:3:length(estimates)]);
        if tshap(1)==1
            ttshap = ['exp         '];
        else
            ttshap = ['ray         '];
        end
        for t=2:length(tshap)
            if(tshap(t)==2)
                ttshap = [ttshap,'ray         '];
            else
               ttshap = [ttshap,'exp         '];
            end
        end
        textstring{1} = (['dist  =  ', num2str(ttshap)]);
        textstring{2} = (['amp   =  ', num2str(tramps)]);
        textstring{3} = (['sig   =  ', num2str(tsigs)]);
        text(50,20,textstring);
        pause(0.05);
        %figure
        BICtemp = BICvalue;
        
        shaps = '';
        for ps=1:length(tshap)
            if tshap(ps)==1
                shaps = [shaps,'E'];
            else
                shaps = [shaps,'R'];
            end
        end
        
        disp(['current minimum: ',shaps,'   BIC=',num2str(BICvalue)]);
    end
    
    % =====     display results: comment/uncomment as desired
    % ====================================================================
    
    
end % of for r
close(h1);


res = fitresults;


end % of function