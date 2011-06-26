% getGaussianPSFsigmaFromData returns the standard deviation of the Gaussian
% approximation of the PSF estimated from point sources in the data.
%
% INPUTS     data     : data structure
%            {frames} : list of frames to include in the estimation. Default: 1.

% Francois Aguet, September 2010


function sigma = getGaussianPSFsigmaFromData(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(img, varargin{:});


% First pass with fixed sigma
pstruct = pointSourceDetection(img, 1.5, 'Mode', 'xyac');
pstruct = fitGaussians2D(img, pstruct.x, pstruct.y, pstruct.A, 1.5*ones(1,length(pstruct.x)), pstruct.c, 'xyasc');

isPSF = pstruct.pval_KS > 0.05 & [pstruct.pval_Ar] > 0.95;
fprintf('PSFs detected: %d\n', sum(isPSF));

svect = pstruct.s(~isnan(pstruct.s) & isPSF);

opts = statset('maxIter', 200);

BIC = zeros(1,3);
sigma = zeros(1,3);
for n = 1:3
    obj = gmdistribution.fit(svect', n, 'Options', opts);
    BIC(n) = obj.BIC;
    sigma(n) = obj.mu(obj.PComponents==max(obj.PComponents));
end

sigma = sigma(BIC==min(BIC));

if strcmpi(ip.Results.Display, 'on')
    ds = 0.2;
    si = -1:ds:10;
    ni = hist(svect, si);
    ni = ni/sum(ni*ds);
    figure;
    h = bar(si,ni);
    set(h, 'BarWidth', 1);
    
    hold on;
    plot(si, pdf(obj,si'), 'r'); 
end
