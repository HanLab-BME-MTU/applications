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


function [prmVect, residual, estimatesSigma, BIC] = fitNWeibull(t, data, prmVect, estVect, mode, display)

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
    'TolX', 1e-12, ...
    'Tolfun', 1e-12);


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


% % degrees of freedom = number of data points minus number of free fit parameters
% numFreeFitP = sum(estVect);
% degFreedom  = numel([data{:}]) - numFreeFitP;
% % chi square = sum of residual divided by degrees of freedom
% chiSquare   = nansum(residual.^2)/degFreedom;
% % cofactor matrix Q; since inverse on JJ isn't possible because of the
% % zeros at positions of fixed parameters, perform the inverse operation on
% % a condensed version of JJ with no zeros
% JJ          = full(jacobian)'*full(jacobian);
% JJdefpos    = find(JJ~=0);
% if length(JJdefpos)==numFreeFitP^2
%     JJsmall     = zeros(numFreeFitP);
%     JJsmall(:)  = JJ(JJdefpos);
%     Qsmall      = inv(JJsmall);
%     Q           = zeros(length(estVect));
%     Q(JJdefpos)    = Qsmall(:);
%     
%     % standard deviation of parameters uses only diagonal of covariance matrix,
%     % which is cofactor times error
%     % Note: large values outside of the diagonal (i.e. on the same order of
%     % magnitude as the diagonal) indicate interdependence of parameters
%     estimatesSigma = sqrt(chiSquare*diag(Q))';
% else
%     estimatesSigma = [];
% end
estimatesSigma = [];


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