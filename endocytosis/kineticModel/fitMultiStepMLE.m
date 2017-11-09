%[k, y, BIC] = fitMultiStepMLE(x, f, varargin) implements a multi-step process 
% of the form 
%
%    k1     k2      kn
% S0 --> S1 --> ... --> Sn
%
% INPUTS
%         x : histogram bin boundaries
%      data : histogram counts
%
% OPTIONS (Specifier, 'value')
%     'Display' : true|{false} plots the fitting result
%       'Noise' : 'Gaussian'|{'Poisson'} selects the noise model
%
% OUTPUTS
%         k : estimated rates
%    k_pstd : s.d. of rates, obtained by error propagation
%       res : structure with fields:
%             .Lmax : maximized log-likelihood value and s.d.
%             .BIC : Bayesian Information Criterion and s.d.
%             .AICc : corrected Akaike Information Criterion and s.d.
%             .model : fitted model. If 'CalculateModelUncertainty' was 'true'
%                      the second column contains the model s.d.
%             .C : covariance matrix of the parameters
%
% Notes: FFT-based implementation

% Francois Aguet, 2013-2014

function [k, k_pstd, res, BIC] = fitMultiStepMLE(x, data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addOptional('StepRange', 1:5);
ip.addParamValue('Noise', 'Poisson', @(x) any(strcmpi(x, {'Gaussian', 'Poisson'})));
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('NumInit', 1, @isposint);
ip.addParamValue('Reference', []);
ip.parse(x, data, varargin{:});
stepRange = ip.Results.StepRange;

K = numel(stepRange);
p = cell(1,K);

% single-step (exp) fit
mu = sum(data.*x)/sum(data);
k0 = 1/mu;
[p{1}, res(1)] = fitHistMLE(x, data, @stepModelFFT, k0, 'Noise', ip.Results.Noise);

nIter = ip.Results.NumInit;
for i = 2:K
    % perform fit multiple times with random intializations
    pv = cell(1,nIter);
    resv = cell(1,nIter);
    for j = 1:nIter
        % update initial values
        k = p{i-1}(:,1);
        k0 = [k; (rand*0.8+0.1)*max(k)];
        %k0 = rand(1,i)*0.05;
        [pv{j}, resv{j}] = fitHistMLE(x, data, @stepModelFFT, k0, 'Noise', ip.Results.Noise);
    end
    % select best fit
    resv = vertcat(resv{:});
    BIC = vertcat(resv.BIC);
    [~,idx] = min(BIC(:,1));
    p{i} = pv{idx};
    res(i) = resv(idx);
end

% best fit
BIC = vertcat(res.BIC);

% pairwise t-test between successive BIC values
hval = zeros(1,K-1);
pval = zeros(1,K-1);
for i = 1:K-1
    [hval(i), pval(i)] = ttest2Parameters(BIC(i:i+1,1), BIC(i:i+1,2), numel(data)*[1 1]);
end

% last significant difference:
% [~,i] = min(BIC(:,1));
%idx = find(hval==1, 1, 'last')+1;
idx = find(hval==0, 1, 'first'); % in case of numerical errors with more complex models
if isempty(idx)
    idx = K;
end

% Output:
k = p{idx}(:,1);
k_pstd = p{idx}(:,2);
res = res(idx);

% chi2 test
% s = sum((data(:)-res.model(:,1)).^2./res.model(:,1));
% pchi2 = chi2cdf(s,numel(data)-1);

if ip.Results.Display
    ha = setupFigure(1,3, 'XSpace', [2 2 0.5], 'YSpace', [1.5 0.5 1]);
    plot(ha(1), x, data, 'k', 'LineWidth', 1);
    plot(ha(1), x, res.model(:,1), 'r', 'LineWidth', 1);
    xlabel(ha(1), 'Lifetime (s)');
    ylabel(ha(1), '# Events');
    set(ha(1), 'XLim', [0 240], 'XTick', 0:40:240);
    if ~isempty(ip.Results.Reference)
        plot(ha(1), x, ip.Results.Reference, 'b--');
    end
    
    % rates
    he = errorbar(ha(2), k, k_pstd, 'LineStyle', 'none',...
        'Color', 0.6*[1 1 1], 'LineWidth', 2);
    setErrorbarStyle(he, 0);
    plot(ha(2), 1:idx, k, 'k.', 'MarkerSize', 12);
    set(ha(2), 'XLim', [0.5 idx+0.5], 'XTick', 1:3);
    YLim = get(ha(2), 'YLim');
    set(ha(2), 'YLim', [0 YLim(2)]);
    xlabel(ha(2), 'Step #');
    ylabel(ha(2), 'Rate (s^{-1})');
    
    % BIC
    he = errorbar(ha(3), BIC(:,1), BIC(:,2), 'LineStyle', 'none',...
        'Color', 0.6*[1 1 1], 'LineWidth', 2);
    setErrorbarStyle(he, 0);
    plot(ha(3), BIC(:,1), 'k.', 'MarkerSize', 12);
    set(ha(3), 'XLim', [0.5 K+0.5], 'XTick', 1:K);
    % indicate significant difference
    YLim = get(ha(3), 'YLim');
    i = idx-1;
    plot(ha(3), i+[0.05 0.05 0.95 0.95],...
        max(BIC(i:i+1,1)+BIC(i:i+1,2))+diff(YLim)*[0.04 0.08 0.08 0.04],...
        'k', 'LineWidth', 1);
    text(i+0.5, BIC(i,1)+BIC(i,2)+diff(YLim)*0.1, ['p < 10^{' num2str(ceil(log10(pval(i)))) '}'], 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center', 'Parent', ha(3));
    % indicate n.s. differences
    for i = idx:K-1
        if hval(i)==0
            plot(ha(3), i+[0.05 0.05 0.95 0.95],...
                max(BIC(i:i+1,1)+BIC(i:i+1,2))+diff(YLim)*[0.04 0.08 0.08 0.04],...
                'k', 'LineWidth', 1);
            text(i+0.5, BIC(i,1)+BIC(i,2)+diff(YLim)*0.1, 'n.s.', 'VerticalAlignment', 'bottom',...
                'HorizontalAlignment', 'center', 'Parent', ha(3));
        end
    end
    xlabel(ha(3), 'Model #');
    ylabel(ha(3), 'BIC');
    
    %K = corrMatFromCov(res(idx).C);
    %plotCorrelationMatrix(K, 'TickLabels', arrayfun(@(i) ['k_' num2str(i)], 1:i, 'unif', 0));
end
