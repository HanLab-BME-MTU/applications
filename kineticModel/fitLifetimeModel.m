%res = fitLifetimeModel(lftData, varargin)
%
% Inputs:
%         lftData : structure returned by 'runLifetimeAnalysis'
%
% Options:
%          'Mode' : 'PDF | {'CDF'} fit to the histogram or empirical distribution
%          'NumP' : Number of populations to test/fit. Default: 1:3
%  'ConstrainBIC' : Smallest BIC with probability > AlphaBIC than next candidate is chosen
%      'AlphaBIC' : Threshold for BIC selection; default: 0.95
%       'PlotAll' : Display correlation matrix and BIC values
%       'PlotCDF' : Display CDF and fitted model
%
% Output:
%             res : output structure with fields:
%                .k            : vector of rates from the model
%                .k_std        : standard deviation (error propagated) of the rates
%                .corr         : correlation matrix for the rates
%                .BIC          : BIC for the populations tested
%                .pPercentiles : percentiles of each subpopulations in the optimal model
%                .pA           : constributions of each subpopulation in the optimal model
%                .pMean        : means of each subpopulation in the optimal model

% Francois Aguet (last modified 01/23/2012)

function res = fitLifetimeModel(lftData, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftData');
ip.addParamValue('Mode', 'CDF', @(x) any(strcmpi(x, {'PDF', 'CDF'})));
ip.addParamValue('NumP', 3, @(x) all(ismember(x, 1:4)));
ip.addParamValue('ConstrainBIC', true, @islogical);
ip.addParamValue('AlphaBIC', 0.95);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('PlotCDF', false, @islogical);
ip.addParamValue('Display', true, @islogical);
ip.addParamValue('Verbose', false, @islogical);
ip.addParamValue('JackKnife', false, @islogical);
% ip.addParamValue('ShowInset', false, @islogical);
ip.addParamValue('EndIdx', find(lftData.meanHist~=0, 1, 'last'));
ip.parse(lftData, varargin{:});
dBIC = 2*log(ip.Results.AlphaBIC/(1-ip.Results.AlphaBIC));


% cutoff at last non-zero data point
endIdx = ip.Results.EndIdx;

t = lftData.t(1:endIdx);
dt = t(2)-t(1);

lftECDF = cumsum(lftData.meanHist)*dt;
lftECDF = lftECDF(1:endIdx);
a = lftECDF(1);
lftECDF = (lftECDF-a)/(1-a);

lftHist = lftData.meanHist(1:endIdx);
lftHist = lftHist/sum(lftHist)/dt;

%===================================================
% Weights for weighted least-squares fit
%===================================================
switch ip.Results.Mode
    case 'PDF'
        W = std(vertcat(lftData.lftHist{:}), [], 1);
        %W = 1/norminv(0.75) * mad(vertcat(lftData.lftHist{:}), 1, 1);
        W = W(1:endIdx);
        W = 1./W;
        W = W/max(W);
    case 'CDF'

        lftECDF_all = cumsum(vertcat(lftData.lftHist{:}), 2);
        lftECDF_all = lftECDF_all(:,1:endIdx);
        A = lftECDF_all(:,1);
        for i = 1:numel(A)
            lftECDF_all(i,:) = (lftECDF_all(i,:)-A(i))/(1-A(i));
        end
        W = std(lftECDF_all, [], 1);
        %W = (1/norminv(0.75) * mad(lftECDF_all, 1, 1)).^2;
        
        W(W==0) = min(W(W~=0));
        W = 1./W; % weight: 1/sigma^2 -> 1/sigma for lsqnonlin
        W = W/max(W);
end
if ip.Results.PlotAll
    fset = loadFigureSettings();

    figure;
    hold on;
    plot(t, W, 'r-', 'LineWidth', 1.5);
    %set(gca, 'LineWidth', 1.5, 'Layer', 'top', fset.sfont{:}, 'XLim', [t(1) t(end)]);
    set(gca, 'LineWidth', 1.5, 'Layer', 'top', fset.sfont{:});
    axis([t(1) t(end) 0 1]);
    xlabel('t', fset.lfont{:});
    ylabel('W(t)', fset.lfont{:});
end







opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', 1e5, ...
    'Display', 'off', ...
    'TolX', 1e-12, ...
    'Tolfun', 1e-12);

dti = dt/10;
t_fine = 0:dti:t(end);
n = numel(lftHist);

maxp = max(ip.Results.NumP);
res.k = cell(1,maxp);
res.k_pstd = cell(1,maxp);
res.corr = cell(1,maxp);
res.BIC = zeros(1,maxp);

for i = ip.Results.NumP
    
    % # states
    ns = i*2;
    
    % Intializations & bounds
    switch i
        case 1
            k0 = 0.01;
        case 2
            k0 = [0.2 0.2 0.05];
        case 3
            %k0 = [0.2 0.2 0.05 0.01 0.02];
            k0 = [0.6    0.25    0.1    0.1    0.03];
            %k0 = 0.1 * ones(1,ns-1);
        case 4
            %k0 = [0.2 0.2 0.05 0.01 0.02 0.01 0.02];
            k0 = 0.02 * ones(1,ns-1);
    end
    lb = zeros(1,ns-1);
    ub = Inf(1,ns-1);
    
    switch ip.Results.Mode
        case 'PDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costPDF, k0, lb, ub, opts, t, lftHist, i, W);
%             k = fminsearch(@costPDFmedian, k0, opts, t, lftHist, i);
%           
%             %opts.('Algorithm') = 'levenberg-marquardt';
%             %opts.('LargeScale') = 'off';
%             %[k, fval, exitflag, output, grad, hessian] = fminunc(@(x) costPDFmedian(x, t, lftHist, i), k0, opts);
%         
%             %opts.('Algorithm') = 'sqp';
%             %[k,fval,exitflag,output,lambda,grad] = fmincon(@(x) costPDFmedian(x, t, lftHist, i), k0, [], [], [], [], lb, ub, [], opts)
%             
%             residual = costPDF(k, t, lftHist, i);
%             resnorm = sum(residual.^2);
%             
%             h = 2*sqrt(eps); % assuming that f''(x) ~ 1
%             nk = numel(k);
%             J = zeros(numel(t), nk);
%             for ki = 1:nk
%                 kp = k;
%                 kp(ki) = kp(ki)+h;
%                 J(:,ki) = (computePDF(kp, t, i) - computePDF(k, t, i))/h;
%             end
        case 'CDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costCDF, k0, lb, ub, opts, t, lftECDF, i, W);
            %k = fminsearch(@costCDFmedian, k0, opts, t, lftECDF, i);
            %residual = costCDF(k, t, lftECDF, i);
            %resnorm = sum(residual.^2);
            
            %h = 2*sqrt(eps); % assuming that f''(x) ~ 1
            %nk = numel(k);
            %J = zeros(numel(t), nk);
            %for ki = 1:nk
            %    kp = k;
            %    kp(ki) = kp(ki)+h;
            %    J(:,ki) = (computeCDF(kp, t, i) - computeCDF(k, t, i))/h;
            %end
    end
    %k
    % BIC and correlation matrix
    J = full(J);
    C = resnorm/(n-numel(k)-1)*inv(J'*J);
    k_pstd = sqrt(diag(C))';
    K = corrFromC(C)';
    
    res.k{i} = k;
    res.k_pstd{i} = k_pstd;
    res.corr{i} = K;
    res.BIC(i) = n*log(resnorm/n) + numel(k)*log(n);
end

% Only significant differences in BIC (alpha = 0.05 default) are considered
if ip.Results.ConstrainBIC && numel(res.BIC)>1
    sortBIC = sort(res.BIC);
    minIdx = find(res.BIC==min(sortBIC(diff(sortBIC) > dBIC)));
else
    minIdx = find(res.BIC==min(res.BIC));
end
np = minIdx;

if ip.Results.JackKnife
    % bootstrap optimal p
    N = numel(lftData.lftHist);
    M = vertcat(lftData.lftHist{:});
    M = M(:,1:endIdx);
    k_jk = cell(1,N);
    parfor i = 1:N
        jkMean = mean(M(setdiff(1:N,i),:), 1);
        switch ip.Results.Mode
            case 'PDF'
                k_jk{i} = lsqnonlin(@costPDF, k0, lb, ub, opts, t, jkMean, np);
            case 'CDF'
                jkECDF = cumsum(jkMean)*dt;
                k_jk{i} = lsqnonlin(@costCDF, k0, lb, ub, opts, t, jkECDF, np);
        end
    end
    k_jk = vertcat(k_jk{:});
    k_std = std(k_jk, [], 1) / sqrt(N);
else
    k_std = res.k_pstd{np};
end

% Intializations & bounds
ns = np*2;
S0 = [1 zeros(1,ns-1)];

sol = ode45(@(t,y) getStateMatrix(np, res.k{np})*y, [0 t(end)], S0);
Y = deval(sol, t_fine);
% A = getStateMatrix(np, res.k{np});
% ev = eig(A)

if ip.Results.Verbose
   for i = 1:numel(res.k{np})
       fprintf('k%d = %.4f, ', i, res.k{np}(i));
   end
   fprintf('\b\b\n');
end

% normalize subpopulation PDFs, weigh by output
popPDF = zeros(np,numel(t_fine));
for i = 1:np
    p = Y(2*i-1,:);
    popPDF(i,:) = p/sum(p)/dti * Y(2*i,end);
end
pdf = sum(popPDF, 1);

popCDF = Y(2:2:end,:);
cdf = sum(popCDF,1);
a = interp1(t_fine, cdf, t(1));

% Compute population percentiles, mean, and contribution
pECDF = arrayfun(@(i) Y(i,:)/Y(i,end), 2:2:ns, 'UniformOutput', false);
for i = 1:np
    [u, uidx] = unique(pECDF{i});
    res.pPercentiles{i} = interp1(u, t_fine(uidx), [0.05 0.25 0.5 0.75 0.95]);
end
res.pA = Y(2:2:end,end)';
res.pMean = arrayfun(@(i) sum(t_fine.*popPDF(i,:)*dti) / sum(popPDF(i,:)*dti), 1:np);

res.np = np;
res.a = a;
res.t = t;
res.t_fine = t_fine;
res.pdf = pdf;
res.popPDF = popPDF;
res.cdf = cdf;
res.popCDF = popCDF;
res.k_std = k_std;
res.lftHist = lftHist;
res.lftECDF = lftECDF;

if ip.Results.Display
    plotLifetimeModel(res, 'PlotAll', ip.Results.PlotAll, 'PlotCDF', ip.Results.PlotCDF);
end



function pdf = computePDF(kVect, t, M)

% Intializations & bounds
S0 = [1 zeros(1,2*M-1)];

% interpolate result over full time vector
dt = t(2)-t(1);
t_full = 0:dt:t(end);

sol = ode45(@(t,y) getStateMatrix(M,kVect)*y, [0 t(end)], S0);
Y = deval(sol, t_full); %[Y, dY]

% normalize subpopulation PDFs, weigh by output
popPDF = zeros(M,numel(t_full));
for k = 1:M
    p = Y(2*k-1,:);
    popPDF(k,:) = p/sum(p)/dt * Y(2*k,end);
end
pdf = sum(popPDF, 1);

%normalize pdf to 1 over 't'
pdf = interp1(t_full, pdf, t);
n = sum(pdf)*dt;
pdf = pdf/n;

function v = costPDF(kVect, t, lftHist, M, W)
if nargin<5
    W = ones(size(lftHist));
end
v = computePDF(kVect, t, M) - lftHist;
v = v.*W;

function v = costPDFmedian(kVect, t, lftHist, M)
v = costPDF(kVect, t, lftHist, M);
v = median(v.^2);



function cdf = computeCDF(kVect, t, M)
S0 = [1 zeros(1,2*M-1)];

sol = ode45(@(t,y) getStateMatrix(M,kVect)*y, [0 t(end)], S0);
Y = deval(sol, t);
cdf = sum(Y(2:2:end,:),1);

% normalize to [0..1]
T = cdf(1);
cdf = (cdf-T)/(1-T);

function v = costCDF(kVect, t, lftECDF, M, W)
if nargin<5
    W = ones(size(lftECDF));
end
v = computeCDF(kVect, t, M) - lftECDF;
v = v.*W;


function v = costCDFmedian(kVect, t, lftECDF, M)
v = costCDF(kVect, t, lftECDF, M);
v = median(v.^2);



function K = corrFromC(C)
n = size(C,1);
K = zeros(n,n);

idx = pcombs(1:n);
i = idx(:,1);
j = idx(:,2);
ij = i+n*(j-1);
ii = i+n*(i-1);
jj = j+n*(j-1);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));
