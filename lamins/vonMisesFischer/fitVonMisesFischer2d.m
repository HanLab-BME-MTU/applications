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
% vmf = vonMisesFischer2d(theta);
% vmfres = @(x) data - x(2:2:end-2)*vmf(x(1:2:end-2)',x(end-1)) - x(end);

vmfres = vonMisesFischer2dMult(theta,data);

r = zeros(1,maxN);

init = zeros(1,maxN*2+2);
init(1:2:end-2) = rand(1,maxN)*pi-pi/2;
% init(2:2:end-2) = rand(1,maxN);
init(end-1) = length(data);
init(end) = mean(data);

[maxv,maxres] = max(data);
init(1) = theta(maxres);
init(2) = (maxv-init(end)) / (1./2/pi/besseli(0,init(end-1))*exp(init(end-1)));
% init(end) = mean(data - vmfres(init([1:2 end-1 end])));

muMax = 0;
if(nargin > 3)
    muMax = min(length(mu)*2,length(init)-2);
    init(1:2:muMax) = mu(1:muMax/2);
end

% options = optimoptions(@lsqnonlin,'PlotFcns',@optimplotfval,'TolX',1e-8,'TolFun',1e-8);

options = optimoptions(@lsqnonlin,'TolX',1e-8,'TolFun',1e-8,'Display','off','Jacobian','on');

lb = zeros(1,(maxN+1)*2);
ub = Inf(1,(maxN+1)*2);
lb(2:2:end-2) = 0.1;
lb(1:2:end-2) = -pi/2;
% lb(end-1) = 0.3;
lb(end) = 0;
ub(1:2:end-2) = pi/2;

numDegFree_last = 0;
resnorm_last = 0;
bic_last = 0;

for n=1:maxN
    nparams = (n+1)*2;
    % aT, trial solution
    % resnormT, trial residual
    [aT,resnormT,residual] = lsqnonlin(vmfres,[init(1:nparams-2) init(end-1:end)],[lb(1:nparams-2) lb(end-1:end)],[ub(1:nparams-2) ub(end-1:end)],options);
    fits{n} = aT;
%     r(n) = sum(vmfres(aT).^2);
    r(n) = resnormT;
%     if(r(n) < tol)
%         break;
%     end

    % 
    numDegFreeT = length(data)-length(aT);
    bicT = length(data)*log(resnormT/length(data))+log(length(data))*length(aT);
    bic(n) = bicT;
    
    if(n > 1)
            %get test statistic, which is F-distributed
%             testStat = (sum(residualsT.^2)/numDegFreeT)/...
%                 (sum(residuals.^2)/numDegFree);
            testStat = (resnormT/numDegFreeT)/...
                (resnorm_last/numDegFree_last);
            
            %get p-value of test statistic
            pValue = fcdf(testStat,numDegFreeT,numDegFree_last);
            
            dBic = bicT-min(bic(1:end-1));
            evidenceRatio = exp(0.5*(bicT-bic_last));
            
            
            %compare p-value to alpha
            %1-sided F-test: H0: F=1, H1: F<1
%             pValue
%            if dBic < 10
            if bicT == min(bic)
%             if pValue < alpha %if p-value is smaller, accept this fit
%                 fit = 1; %and attempt another one with an additional kernel
                numDegFree_last = numDegFreeT;
                resnorm_last = resnormT;
                bic_last = bicT;
                a = aT;
%                 init(1:length(aT)-2) = aT(1:end-2);
% %                 
%                 [maxv,maxres] = max(residual);
%                 init(length(aT)-1) = theta(maxres);
% %                 init(length(aT)) = aT(end-2);
% 
%                 init(length(aT)) = (maxv) / (1./2/pi/besseli(0,aT(end-1))*exp(aT(end-1)));
% %                 
% % %                 init(length(aT)-1:2:muMax) = mu(length(aT)/2:muMax/2);
% %                 
%                 init(length(aT)+1:length(aT)+2) = aT(end-1:end);
%                 init(end-1:end) = aT(end-1:end);
                
            else %if p-value is larger, do not accept this fit and exit
%                 fit = 0;
%                 break;
            
            end
%             if(dBic > 20)
%                 break;
%             end
    else
        bic_last = bicT;
        numDegFree_last = numDegFreeT;
        resnorm_last = resnormT;
        a = aT;
        init(1:length(aT)-2) = aT(1:end-2);
%         
%         [~,maxres] = max(residual);
%         init(length(aT)-1) = theta(maxres);
%         init(length(aT)) = aT(end-2);
%                 [maxv,maxres] = max(residual);
%                 init(length(aT)-1) = theta(maxres);
% %                 init(length(aT)) = aT(end-2);
% 
%                 init(length(aT)) = (maxv) / (1./2/pi/besseli(0,aT(end-1))*exp(aT(end-1)));
% %                 
% % %                 init(length(aT)-1:2:muMax) = mu(length(aT)/2:muMax/2);
% %                 
%                 init(length(aT)+1:length(aT)+2) = aT(end-1:end);
% %         
% % %         init(length(aT)-1:2:muMax) = mu(length(aT)/2:muMax/2);
%         init(length(aT)+1:length(aT)+2) = aT(end-1:end);
%         init(end-1:end) = aT(end-1:end);
    end
        init(1:length(aT)-2) = aT(1:end-2);
        [maxv,maxres] = max(residual);
        % next location
        init(length(aT)-1) = theta(maxres);
        % next scaling factor
        init(length(aT)) = (maxv) / (1./2/pi/besseli(0,aT(end-1))*exp(aT(end-1)));
%         init(length(aT)+1:length(aT)+2) = aT(end-1:end);
%         init(length(aT)+1:length(aT)+2) = aT(end-1:end);
        % carry over kappa and offset
        init(end-1:end) = aT(end-1:end);
end

if(nargout > 2)
    r = r(1:n);
    bic = bic(1:n);
end

% figure;
% plot(r);

% figure;
% plot(theta,data);
% hold on;
% % plot(theta,data - vmfres(a));
% plotVMFfit(a);
% hold off;
end

