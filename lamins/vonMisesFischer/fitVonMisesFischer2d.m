function [ a , r, bic, fits] = fitVonMisesFischer2d( data, alpha, maxN, mu)
%fitVonMisesFischer2d Fits multiple vonMisesFischer distributions to data

if(nargin < 2)
    alpha = 0.05;
end

if(nargin < 3)
    maxN = length(data)/2-2;
end

N = length(data);
theta = -pi/2:pi/N:pi/2-pi/N;

vmfres = vonMisesFischer2dMult(theta,data);

stop.TolResnorm = eps;
stop.TolX = eps;
stop.TolFun = eps;
stop.Fcn = @(x,optimValues,state) strcmp(state,'iter') && optimValues.resnorm < stop.TolResnorm;
options = optimoptions(@lsqnonlin,'TolX',eps,'TolFun',eps,'Display','off','Jacobian','on','DerivativeCheck','off','OutputFcn',stop.Fcn);
% besseli = @(a,b) 1;

init = zeros(1,maxN*2+2);
init(end-1) = length(data);
init(end) = mean(data);

[maxv,maxres] = max(data);
init(1) = theta(maxres);
init(2) = (maxv-init(end)) / (1./2/pi/besseli(0,init(end-1))*exp(init(end-1)));

% muMax = 0;
% if(nargin > 3)
%     muMax = min(length(mu)*2,length(init)-2);
%     init(1:2:muMax) = mu(1:muMax/2);
% end


% lower bound
lb = zeros(1,(maxN+1)*2);
% lb(2:2:end-2) = 0.1;
lb(1:2:end-2) = -pi/2;
lb(end) = 0;

% upper bound
ub = Inf(1,(maxN+1)*2);
ub(1:2:end-2) = pi/2;

numDegFree_last = 0;
resnorm_last = 0;
bic_last = 0;
minBic = Inf;

numDegFreeT = 0;
bicT = 0;

accept = true;
firstTime = true;

use.ftest = true;
use.minbic = false;
use.dbic = false;
breakOnReject = false;

if(nargout > 1)
    % debugging
    r = zeros(1,maxN);
    bic = zeros(1,maxN);
    fits = cell(1,maxN);
end

for n=1:maxN
    nparams = (n+1)*2;
    % aT, trial solution
    % resnormT, trial residual
    % resnormT == sum(vmfres(aT).^2
    s = [1:nparams-2 length(init)-1 length(init)];
    [aT,resnormT,residual,exitflag,output,lambda,jacobian] = lsqnonlin(vmfres,init(s),lb(s),ub(s),options);
    

    if(firstTime)
%         accept = true;
        firstTime = false;
        numDegFreeT = N - nparams;
    elseif(use.ftest)
        %compare p-value to alpha
        %1-sided F-test: H0: F=1, H1: F<1
        %get test statistic, which is F-distributed
        numDegFreeT = N - nparams;
        testStat = (resnormT/numDegFreeT)/...
            (resnorm_last/numDegFree_last);
        %get p-value of test statistic
        pValue = fcdf(testStat,numDegFreeT,numDegFree_last);
        accept = pValue < alpha;
    elseif(use.minbic || use.dbic)
        logN = log(N);
        bicT = N*(log(resnormT) - logN) + logN*nparams;
        accept = bicT < bic_last;
        if(use.dbic)
            dBic = bicT-minBic;
            if(accept)
                minBic = bicT;
            end
    %           evidenceRatio = exp(0.5*(bicT-bic_last));
            accept = dBic < alpha;
        end
    end
    
    if(nargout > 1)
        % debugging
        r(n) = resnormT;
        bic(n) = bicT;
        fits{n} = aT;
    end
    
    if(accept)
        a = aT;
        resnorm_last = resnormT;
        numDegFree_last = numDegFreeT;
        bic_last = bicT;
        if(resnormT < stop.TolResnorm)
            break;
        end
%         if(use.ftest)
%             % adjust tolerance to be below ftest limit
%             nextStat = finv(alpha,numDegFree_last+2,numDegFree_last);
%             options.TolFun = min(nextStat *resnorm_last/numDegFree_last*(numDegFree_last-2) /2,options.TolFun);
%         end
    elseif(breakOnReject)
        break;
    end
   
    % initialize the next parameters at the same location and amplitude
    init(1:length(aT)-2) = aT(1:end-2);
    
    kappa = aT(end-1);
    
    % next location will be centered at the maximum residual
    [maxv,maxres] = max(residual);
    init(length(aT)-1) = theta(maxres);
    
%     maxima = [aT(2:2:end-2)*(1./2/pi/besseli(0,kappa)*exp(kappa)) maxv];
    
    % next scaling factor
%     kappa = length(data);
    init(length(aT)) = (maxv) / (1./2/pi/besseli(0,kappa)*exp(kappa));
%     init(2:2:length(aT)) = maxima / (1./2/pi/besseli(0,kappa)*exp(kappa));

    % carry over kappa and offset
    init(end-1:end) = aT(end-1:end);
    init(end-1) = kappa;
    

    
end

if(nargout > 1)
    % debugging
    r = r(1:n);
    bic = bic(1:n);
    fits = fits(1:n);
end

end

