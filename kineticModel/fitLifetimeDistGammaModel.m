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

function res = fitLifetimeDistGammaModel(lftRes, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftData');
ip.addParamValue('Mode', 'CDF', @(x) any(strcmpi(x, {'PDF', 'CDF'})));
ip.addParamValue('MaxP', 3, @(x) all(ismember(x, 1:4)));
ip.addParamValue('ConstrainBIC', true, @islogical);
ip.addParamValue('AlphaBIC', 0.95);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('PlotCDF', false, @islogical);
ip.addParamValue('Display', true, @islogical);
ip.addParamValue('Verbose', false, @islogical);
ip.addParamValue('JackKnife', false, @islogical);
ip.addParamValue('fYLim', []);
ip.addParamValue('rYLim', []);
ip.addParamValue('HistName', 'lftHist_A');
% ip.addParamValue('ShowInset', false, @islogical);
ip.addParamValue('EndIdx', []);
ip.parse(lftRes, varargin{:});
dBIC = 2*log(ip.Results.AlphaBIC/(1-ip.Results.AlphaBIC));


lftHistMat = lftRes.(ip.Results.HistName);
meanHist =  mean(lftHistMat,1);

% cutoff at last non-zero data point
endIdx = ip.Results.EndIdx;
if isempty(endIdx)
    endIdx = find(meanHist~=0, 1, 'last');
end

t = lftRes.t(1:endIdx);
dt = t(2)-t(1);

lftEDF = cumsum(meanHist)*dt;
lftEDF = lftEDF(1:endIdx);
a = lftEDF(1);
lftEDF = (lftEDF-a)/(1-a);

lftHist = meanHist(1:endIdx);
lftHist = lftHist/sum(lftHist)/dt;

if strcmpi(ip.Results.Mode, 'CDF')
    f = lftEDF;
else
    f = lftHist;
end

maxP = ip.Results.MaxP;
prmVect = cell(1,maxP);
BIC = zeros(1,maxP);
for N = 1:maxP
    [prmVect{N} a(N) BIC(N)] = fitGammaMixture(t, f, 'kna', 'N', N, 'FitMode', ip.Results.Mode);
end

N0 = find(BIC==min(BIC), 1, 'first');
prmVect = prmVect{N0};
BIC = BIC(N0);
a = a(N0);

% means
mu = prmVect(2:3:end)./prmVect(1:3:end);
[~,idx] = sort(mu);

% sort parameter vector: increasing means
I = reshape(1:3*N0, [3 N0]);
idx = reshape(I(:,idx), [1 3*N0]);
prmVect = prmVect(idx);

dti = dt/10;
t = 0:dti:t(end);
res.t = t;
[res.CDF res.popCDF] = gammaMixture(t, prmVect, 'Mode', 'CDF');
[res.PDF res.popPDF] = gammaMixture(t, prmVect, 'Mode', 'PDF');
res.BIC = BIC;
res.N = N0;
res.a = a;
res.k = prmVect(1:3:end);
res.n = prmVect(2:3:end);
res.A = prmVect(3:3:end);
res.percentiles = getCDFPercentiles(res.t, res.popCDF, [0.025 0.5 0.975]);
res.FitMode = ip.Results.Mode;
res.ModelType = 'Gamma';
