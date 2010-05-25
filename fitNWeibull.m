% [prmVect] = fitNWeibull(t, data, prmVect, estVect, mode, display)
% 
% 
% Inputs: 
%           t       : time vector
%           data    : measured histogram
%           prmVect : parameter vector of size 3*N+1
%           estVect : binary vector indicating which parameters are estimated
%           mode    : type of distribution used for the fit: PDF, CDF, or SVF
%           display : '1' or 'display' to visualize fitting and results
%
% Structure of the parameter vector: [offset A(1) lambda(1) k(1) ... A(N) lambda(N) k(N)]
%
% If 't' and 'data' are cell arrays, will fit to all data sets simultaneously.

% Francois Aguet, March 2010


function [prmVect, residual, prmSigma, BIC] = fitNWeibull(t, data, prmVect, estVect, mode, display)

if nargin<6
    display = 0;
    h = [];
end

if ~iscell(t) && ~iscell(data)
    t = {t};
    data = {data};
end
estIdx = (estVect==1);

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e6, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);


if (display~=0)
    h = zeros(1,length(data));
    figure;    
    for n = 1:length(data)
        plot(t{n}, data{n}, 'k.');
        hold on;
        h(n) = plot(t{n}, nWeibull(t{n}, prmVect, mode), 'r');
    end
    axis([0 max([t{:}]) 0 1.1*max([data{:}])]);
    xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 16);
    ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 16);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);

end

% run optimization
[p, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(@nWeibullCost, prmVect(estIdx), [], [], opts, t, data, prmVect, estIdx, mode, h);
prmVect(estIdx) = p;
prmVect(2:end) = abs(prmVect(2:end));


%===========================================================
% BIC/Schwarz criterion (assumption: normal errors)
%===========================================================
n = length([data{:}]); % data points
deg = sum(estIdx); % degrees of freedom
BIC = n*log(resnorm/n) + deg*log(n);


% estimated error variance: RSS/(n-p). n: data points, p: degrees of freedom.
sigma2 = resnorm / (numel([data{:}]) - sum(estVect));

% parameter variance is given by the diagonal elements of the variance-covariance matrix
jacobian = full(jacobian);
prmSigma = sqrt(sigma2 * diag(inv(jacobian'*jacobian)));


% plot subpopulations
if (display~=0)
    np = (length(prmVect)-1)/3;
    colorV = ones(np,3);
    cstep = 360/np;
    colorV(:,1) = (0:cstep:360-cstep)/360;
    ti = unique([t{:}]);
    set(gca, 'ColorOrder', hsv2rgb(colorV));
    [w W] = nWeibull(ti , prmVect, mode);
    plot(ti, W);
end



function J = nWeibullCost(p, t, data, prmVect, estIdx, mode, h)
prmVect(estIdx) = p;
prmVect(2:end) = abs(prmVect(2:end));

J = cell(1,length(data));
for n = 1:length(data)
    w = nWeibull(t{n}, prmVect, mode);
    J{n} = data{n} - w;
    if ~isempty(h)
        set(h(n), 'XData', t{n});
        set(h(n), 'YData', w);
        drawnow
    end
end
J = [J{:}];