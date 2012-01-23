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
ip.addParamValue('Mode', 'PDF', @(x) any(strcmpi(x, {'PDF', 'CDF'})));
ip.addParamValue('NumP', 3, @(x) all(ismember(x, 1:4)));
ip.addParamValue('ConstrainBIC', true, @islogical);
ip.addParamValue('AlphaBIC', 0.95);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('PlotCDF', false, @islogical);
ip.addParamValue('Display', true, @islogical);
ip.addParamValue('JackKnife', true, @islogical);
% ip.addParamValue('ShowInset', false, @islogical);
ip.addParamValue('EndIdx', find(lftData.meanHist~=0, 1, 'last'));
ip.parse(lftData, varargin{:});
dBIC = 2*log(ip.Results.AlphaBIC/(1-ip.Results.AlphaBIC));


% cutoff at last non-zero data point
endIdx = ip.Results.EndIdx;

t = lftData.t(1:endIdx);
dt = t(2)-t(1);

lftHist = lftData.meanHist(1:endIdx);
lftHist = lftHist/sum(lftHist)/dt;
% lftECDF = cumsum(lftHist)*dt;
lftECDF = lftData.meanECDF(1:endIdx);

a = lftECDF(1);
lftECDF = (lftECDF-a)/(1-a);

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', 1e5, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

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
    S0 = [1 zeros(1,ns-1)];
    switch i
        case 1
            k0 = 0.01;
        case 2
            k0 = [0.2 0.2 0.05];
        case 3
            k0 = [0.2 0.2 0.05 0.01 0.02];
            %k0 = [0.5935    0.2466    0.0947    0.1129    0.0331];
            %k0 = 0.01 * ones(1,ns-1);
        case 4
            %k0 = [0.2 0.2 0.05 0.01 0.02 0.01 0.02];
            k0 = 0.02 * ones(1,ns-1);
    end
    lb = zeros(1,ns-1);
    ub = Inf(1,ns-1);
    
    switch ip.Results.Mode
        case 'PDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costPDF, k0, lb, ub, opts, t, lftHist, S0, i);
        case 'CDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costCDF, k0, lb, ub, opts, t, lftECDF, S0, i);
    end
   
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
    for i = 1:N
        jkMean = mean(M(setdiff(1:N,i),:), 1);
        switch ip.Results.Mode
            case 'PDF'
                k_jk{i} = lsqnonlin(@costPDF, k0, lb, ub, opts, t, jkMean, S0, np);
            case 'CDF'
                jkECDF = cumsum(jkMean)*dt;
                k_jk{i} = lsqnonlin(@costCDF, k0, lb, ub, opts, t, jkECDF, S0, np);
        end
    end
    k_jk = vertcat(k_jk{:});
    k_std = std(k_jk, [], 1) / sqrt(N);
else
    k_std = res.k_pstd{np};
end
% k = res.k{np};


% Intializations & bounds
ns = np*2;
S0 = [1 zeros(1,ns-1)];

hf = str2func(['pop' num2str(np) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, res.k{np}), [0 t(end)], S0);

% interpolate result over full time vector
popPDF = zeros(np,numel(t_fine));
Y = interp1(t_ode, Y, t_fine);

% normalize subpopulation PDFs, weigh by output
for i = 1:np
    p = Y(:,2*i-1);
    popPDF(i,:) = p/sum(p)/dti * Y(end,2*i);
end
pdf = sum(popPDF, 1);

popCDF = Y(:,2:2:end)';
cdf = sum(popCDF,1);
a = interp1(t_fine, cdf, t(1));

% Compute population percentiles, mean, and contribution
pECDF = arrayfun(@(i) Y(:,i)/Y(end,i), 2:2:ns, 'UniformOutput', false);
for i = 1:np
    [u, uidx] = unique(pECDF{i});
    res.pPercentiles{i} = interp1(u, t_fine(uidx), [0.05 0.25 0.5 0.75 0.95]);
end
res.pA = Y(end,2:2:end);
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

if ip.Results.Display
    plotLifetimeModel(res, 'PlotAll', ip.Results.PlotAll, 'PlotCDF', ip.Results.PlotCDF);
end



% M: model #
function v = costPDF(kVect, t, lftHist, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);

% interpolate result over full time vector
dt = t(2)-t(1);
t_full = 0:dt:t(end);

popPDF = zeros(M,numel(t_full));
Y = interp1(t_ode, Y, t_full);

% normalize subpopulation PDFs, weigh by output
for k = 1:M
    p = Y(:,2*k-1);
    popPDF(k,:) = p/sum(p)/dt * Y(end,2*k);
end
pdf = sum(popPDF, 1);

%normalize pdf to 1 over 't'
pdf = interp1(t_full, pdf, t);
n = sum(pdf)*dt;
pdf = pdf/n;

v = pdf - lftHist;



function v = costCDF(kVect, t, lftECDF, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);

% interpolate ODE output to input grid
CDF = interp1(t_ode, sum(Y(:,2:2:end),2), t);

% normalize to [0..1]
T = CDF(1);
CDF = (CDF-T)/(1-T);

v = CDF - lftECDF;




% Model:
%       k1
%    S1 --> S2
function dy = pop1Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -k(1)*y(1);
dy(2) = k(1)*y(1);


% Model:
%       k2     k3
%    S1 --> S3 --> S4
% k1 |
%    v
%    S2
function dy = pop2Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -k(3)*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);


% Model:
%       k2     k4     k5
%    S1 --> S3 --> S5 --> S6
% k1 |   k3 |
%    v      v
%    S2     S4
function dy = pop3Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -(k(3)+k(4))*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);
dy(5) = -k(5)*y(5) + k(4)*y(3);
dy(6) = k(5)*y(5);


% Model:
%       k2     k4     k6     k7
%    S1 --> S3 --> S5 --> S7 --> S8
% k1 |   k3 |   k5 |
%    v      v      V
%    S2     S4     S6
function dy = pop4Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -(k(3)+k(4))*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);
dy(5) = -(k(5)+k(6))*y(5) + k(4)*y(3);
dy(6) = k(5)*y(5);
dy(7) = -k(7)*y(7) + k(6)*y(5);
dy(8) = k(7)*y(7);


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
