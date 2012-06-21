% Francois Aguet (last modified 10/30/2011)

function plotTrackFirstFrameIntensities(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Cutoff_f', 4, @isscalar);
ip.addParamValue('Print', false, @islogical);
ip.addParamValue('YLim', [0 100]);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 125 150]);
ip.parse(varargin{:});

cutoff_f = ip.Results.Cutoff_f;

minLft = cutoff_f*data.framerate;
cohortBounds = ip.Results.CohortBounds_s;
cohortBounds(cohortBounds<=minLft) = [];
cohortBounds = [minLft cohortBounds data.movieLength*data.framerate];

mCh = strcmp(data.source, data.channels);

lftData = getLifetimeData(data);

trackLengths_f = lftData.trackLengths(lftData.catIdx==1);

% # cohorts
nc = numel(cohortBounds)-1;

cohortLabels = arrayfun(@(i) [' ' num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-data.framerate) ' s'], 1:nc, 'UniformOutput', false);

kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background

YLim = ip.Results.YLim;


%----------------------------------------------------
% Start buffer, all tracks
%----------------------------------------------------
nb = 5;
pct1 = kLevel*prctile(lftData.sigma_r_Ia(:,1:nb+1,mCh), [50 25 75 2.5 97.5], 1);
pct2 = prctile([lftData.startBuffer_Ia(:,:,mCh) lftData.intMat_Ia(:,1,mCh)], [50 25 75 2.5 97.5], 1);

xlabels = [arrayfun(@(i) [num2str(i*data.framerate) ' s'], -nb:-1, 'UniformOutput', false) {'First frame'}];
figure;
plot([0 4*nb], [0 0], 'k--', 'HandleVisibility', 'off');
boxplot2({pct1, pct2}, ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1,...
    'XLabels', xlabels,...
    'YLabel', 'Intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'North');
set(hl, 'Box', 'off');

if ip.Results.Print 
    print('-depsc2', 'StartBufferIntensities.eps');
end


%----------------------------------------------------
% End buffer, all tracks
%----------------------------------------------------
nt = size(lftData.intMat_Ia,1);
M = zeros(nt,nb+1);
endVal = zeros(nt,1);
for k = 1:nt
    M(k,:) = lftData.sigma_r_Ia(k,trackLengths_f(k)+(0:nb),mCh);
    endVal(k) = lftData.intMat_Ia(k,trackLengths_f(k),mCh);
end
pct1 = kLevel*prctile(M, [50 25 75 2.5 97.5], 1);
pct2 = prctile([endVal lftData.endBuffer_Ia(:,:,mCh)], [50 25 75 2.5 97.5], 1);

xlabels = [{'Last frame'} arrayfun(@(i) ['+' num2str(i*data.framerate) ' s'], 1:5, 'UniformOutput', false)];
figure;
plot([0 4*nb], [0 0], 'k--', 'HandleVisibility', 'off');
boxplot2({pct1, pct2}, ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1,...
    'XLabels', xlabels,...
    'YLabel', 'Intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'North');
set(hl, 'Box', 'off');

if ip.Results.Print 
    print('-depsc2', 'EndBufferIntensities.eps');
end

cPct = zeros(5,nc);
cPctBG = zeros(5,nc);
bPct = zeros(5,nc);
bPctBG = zeros(5,nc);
intRatioPercentiles = zeros(5,nc);
for k = 1:nc
    cIdx = cohortBounds(k)<=trackLengths_f & trackLengths_f<cohortBounds(k+1);
    
    %A0 = arrayfun(@(t) t.A(mCh,1), tracks(cIdx));
    A0 = lftData.intMat_Ia(cIdx,1,mCh);
    
    %maxInts = arrayfun(@(t) max(t.A(mCh,:)), tracks(cIdx));
    maxInts = nanmax(lftData.intMat_Ia(cIdx,:),[],2);
    % percentiles of the max. int. to 1st frame int. ratio
    intRatioPercentiles(:,k) = prctile(maxInts ./ A0, [50 25 75 2.5 97.5]);
    
    bPct(:,k) = prctile(lftData.startBuffer_Ia(cIdx,end,mCh), [50 25 75 2.5 97.5]);
    bPctBG(:,k) = kLevel*prctile(lftData.sigma_r_Ia(cIdx,nb+1,mCh), [50 25 75 2.5 97.5]);
    
    cPct(:,k) = prctile(lftData.intMat_Ia(cIdx,1,mCh), [50 25 75 2.5 97.5]);
    cPctBG(:,k) = kLevel*prctile(lftData.sigma_r_Ia(cIdx,nb+1,mCh), [50 25 75 2.5 97.5]);
end

figure;
boxplot2({cPctBG, cPct}, ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1,...
    'XLabels', cohortLabels,...
    'YLabel', 'Intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'North');
set(hl, 'Box', 'off');
if ip.Results.Print 
    print('-depsc2', '1stFrameIntensitiesPerCohort.eps');
end


figure;
boxplot2({bPctBG, bPct}, ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1,...
    'XLabels', cohortLabels,...
    'YLabel', 'Intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'North');
set(hl, 'Box', 'off');
if ip.Results.Print 
    print('-depsc2', 'LastBufferFrameIntensityPerCohort.eps');
end


% Plot ratio between the maximum and 1st frame intensities >>> relative size of the objects
% figure('Position', [440 378 640 320], 'PaperPositionMode', 'auto');
% barplot2(meanMaxInt, stdMaxInt, 'ErrorbarPosition', 'top', ...
%     'FaceColor', [0.7 0.9 1], 'EdgeColor', [0 0.7 1],...
%     'GroupDistance', 1, 'BarWidth', 1, 'Angle', 45,...
%     'XLabels', cohortLabels,...
%     'XLabel', 'Lifetime cohort', 'YLabel', 'Max. / 1^{st} frame intensity');

figure('Position', [440 378 640 320], 'PaperPositionMode', 'auto');
boxplot2(intRatioPercentiles, ...
    'FaceColor', [0.7 0.9 1], 'EdgeColor', [0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1,...
    'XLabels', cohortLabels,...
    'XLabel', 'Lifetime cohort', 'YLabel', 'Max. / 1^{st} frame intensity', 'YLim', [0 8]);
if ip.Results.Print 
    print('-depsc2', 'MaxToFirstIntensityRatioCohort.eps');
end

