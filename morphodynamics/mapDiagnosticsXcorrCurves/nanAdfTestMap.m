function [pvec, adfMap] = nanAdfTestMap(map, mapName, smParam)
% nanAdfTestMap Perform an augmented Dickey-Fuller (ADF) test to check the
% stationarity of time series of each window and Plot the result.
%
% Output:
%       pvec        - P-values. The null hypothesis is that the time series
%                   is non-stationary. Thus, a small P-value indicates
%                   stationarity.
%       adfMap      - a matlab figure object in which nonstationary windows
%                   are shown with no transparency, and stationary windows
%                   are shown with transparency. 
%
% Jungsik Noh, 2016/10/05


figFlag = 'off';

indNanRows = any(isnan(map'))';
map2 = map(~indNanRows,:);

%map2 = knnimpute(map')';

% adf
tmax = size(map2, 2);
% automatic lag selection
lag0 = floor((tmax-1)^(1/3));
pvec = NaN*ones(size(map2, 1), 1);

for w = 1:size(map2, 1)
    x = map2(w, :);
    %[h, pval, stat] = adftest(x, 'model', 'TS', 'lags', lag0);
    [~, pval, ~] = adftest(x, 'model', 'TS', 'lags', lag0);
    pvec(w) = pval;
end

adfIndex = (pvec < 0.1);           % stationary: 1, nonstationary: 0

% Transform to original map's size
adfIndexOrig = zeros(size(map, 1), 1);
adfIndexOrig(~indNanRows) = adfIndex;

alphaIndex = (1 - 0.8*adfIndexOrig);    % Nonstationary -> alphaIndex=1 (shown)
%
alphaIndex(indNanRows) = 0;

filteredmap = smoothActivityMap(map, 'SmoothParam', smParam, 'UpSample', 1);
adfMap = figure('Visible', figFlag);
figtmp = imagesc(filteredmap);
colorbar;colormap(jet)
figtmp.AlphaData = repmat(alphaIndex, 1,  tmax);
axis xy;xlabel('Time frame');ylabel('Window')
title({mapName, ['total num: ', num2str(size(map2,1)), ' num of nonstationary TS: ', num2str(sum(~adfIndex))] } ) 

end


% ftn:end
