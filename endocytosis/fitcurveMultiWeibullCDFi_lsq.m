function [estimates, residual, estimatesSigma, BICvalue] = fitcurveMultiWeibullCDFi_lsq(t, data, guessVec, fixVec, weightvec)
% fit with inverse of multi-Weibull cumulative distribution function
% the number of weibulls fitted is determined by the number of points in
% the startpoint guess vector (guessVec)
% SYNOPSIS:
% INPUT     t
%           data
%           guessVec
%           fixvec
% OUTPUT    estimates
%           residual
%           estimatesSigma
%           BICvalue        = Bayesian information criterion

% guessVec contains, in that order:
% offset, a1, lambda1, k1, a2, lambda2, k2, ... and so forth

start_point = guessVec; 

% number of Weibulls
num = round((length(guessVec)-1)/3);

fixvector = zeros(3*num+1,1);

if mod(length(guessVec),3)~=1
    error('guessVec has the wrong number of points, it should be 3*n+1');
end

% keep variables fixed? 
% default: no fixing
fixvar = 0;
if nargin>3
    % if fixVec is a vector and not all zeros
    if length(fixVec)>1 & max(fixVec)>0 
        fixvar = 1;
        fixvector = zeros(3*num+1,1);
        fixvector(1:length(fixVec)) = fixVec;
    end
end

% weighting
wei = 0;
if nargin>4
    if length(weightvec)==length(data)
        wei=1;
    end
end
    

options = optimset('Display','off');

%estimates = fminsearch(@expfun, start_point);
[estimates,renorm,residual,exFlag,outP,lam,jacobian]= lsqnonlin(@multipleWeibullCDFun, start_point,[],[],options);

% user-defined function:
% in lsqnonlin: instead of sse, output function value itself, since the 
% squared sum is done implicitly

    function F = multipleWeibullCDFun(params)
        
        % check if there's fixing, and if so reset parameters to the
        % appropriate values of the startvector
        if fixvar == 1
            for f=1:length(fixvector)
                if fixvector(f)==1
                    params(f) = start_point(f);
                end % of if
            end % of for
        end % of if       
        
        % set all parameters except offset to positive
        params(2:length(params))=abs(params(2:length(params)));
        
             
        [FittedCurve]=multiWeibullCDFi(t,params);
           
        % include weighting here if specified
        if wei==1
            F = weightvec.*(FittedCurve - data);
        else
            F = FittedCurve - data;
        end
        
                
        % plot while fitting
        v = axis;
        plot(t,data,'b-');
        hold on
        plot(t,FittedCurve,'r-');
        hold off
        axis(v);
        pause(0.01);     
                
    end


% set final result for estimate to positive
estimates(2:length(estimates))=abs(estimates(2:length(estimates)));

%=========================================================================
%                               fit statistics

% degrees of freedom = number of data points minus number of free fit
% parameters
numFitP     = length(fixvector);
numFreeFitP = length(find(fixvector==0));
numDataP    = sum(isfinite(data));
degFreedom  = numDataP - numFreeFitP;
% chi square = sum of residuals divided by degrees of freedom
chiSquare   = nansum(residual.^2)/degFreedom;
% cofactor matrix Q; since inverse on JJ isn't possible because of the
% zeros at positions of fixed parameters, perform the inverse operation on
% a condensed version of JJ with no zeros
JJ          = full(jacobian)'*full(jacobian);
JJdefpos    = find(JJ~=0);
if length(JJdefpos)==numFreeFitP^2
    JJsmall     = zeros(numFreeFitP);
    JJsmall(:)  = JJ(JJdefpos);           
    Qsmall      = inv( JJsmall ) ;
    Q           = zeros(numFitP);
    Q(JJdefpos)    = Qsmall(:);

% standard deviation of parameters uses only diagonal of covariance matrix,
% which is cofactor times error
% Note: large values outside of the diagonal (i.e. on the same order of
% magnitude as the diagonal) indicate interdependence of parameters
    estimatesSigma = sqrt(chiSquare*diag(Q))';
else
    estimatesSigma = [];
end


%=========================================================================
%                               BIC value
% ( Bayesian information criterion )
    
% n = sample size
BIC_n = numDataP;
% k = number of free parameters
BIC_k = numFreeFitP;
% rss = residual sum of squares
BIC_rss = sum(residual.^2);
    
BICvalue = (BIC_n*log(BIC_rss/BIC_n)) + BIC_k*log(BIC_n);
        



%=========================================================================
% plot final results: individual populations
tplot = [0.5:0.5:max(t)];
hold on
p1params = [0 estimates(2:4)];
p1curve = multiWeibullCDFi(tplot,p1params);
plot(tplot,p1curve,'c-');
if length(estimates)>6
    p2params = [0 estimates(5:7)];
    p2curve = multiWeibullCDFi(tplot,p2params);
    plot(tplot,p2curve,'g-');
    if length(estimates)>9
        p3params = [0 estimates(8:10)];
        p3curve = multiWeibullCDFi(tplot,p3params);
        plot(tplot,p3curve,'m-');
        if length(estimates)>12
            p4params = [0 estimates(11:13)];
            p4curve = multiWeibullCDFi(tplot,p4params);
            plot(tplot,p4curve,'r-');
        end
            
    end
end


end