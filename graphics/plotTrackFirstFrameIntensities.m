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



%%
%----------------------------------------------------
% Start buffer, all tracks
%----------------------------------------------------
nb = 5;
nIa = size(lftData.intMat_Ia,1);
bufferMean = mean(lftData.startBuffer_Ia(:,:,mCh),1);
bufferStd = std(lftData.startBuffer_Ia(:,:,mCh),[],1) / sqrt(nIa);
bufferBGMean = mean(lftData.sigma_r_Ia(:,1:5,mCh),1);
bufferBGStd = std(lftData.sigma_r_Ia(:,1:5,mCh),[],1) / sqrt(nIa);

% intensity in the first frame
A0 = lftData.intMat_Ia(:,1,mCh);
bufferMean(nb+1) = mean(A0);
bufferStd(nb+1) = std(A0)/sqrt(numel(A0));
c0 = lftData.sigma_r_Ia(:,1,mCh);
bufferBGMean(nb+1) = mean(c0);
bufferBGStd(nb+1) = std(c0)/sqrt(numel(c0));
  
xlabels = [arrayfun(@(i) [num2str(i*data.framerate) ' s'], -5:-1, 'UniformOutput', false) {'First frame'}];
figure;
barplot2([kLevel*bufferBGMean; bufferMean]', [kLevel*bufferBGStd; bufferStd]', 'ErrorbarPosition', 'both', ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1, 'YLim', YLim,...
    'XLabels', xlabels,...
    'YLabel', 'Mean intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'NorthOutside');
set(hl, 'Box', 'off');

if ip.Results.Print 
    print('-depsc2', 'StartBufferIntensities.eps');
end

%%
%----------------------------------------------------
% End buffer, all tracks
%----------------------------------------------------
% bufferMean = zeros(5,1);
% bufferStd = zeros(5,1);
% bufferBGMean = zeros(5,1);
% bufferBGStd = zeros(5,1);
% for k = 1:5
%     A0 = arrayfun(@(t) t.endBuffer.A(mCh,k), tracks);
%     bufferMean(k+1) = mean(A0);
%     bufferStd(k+1) = std(A0)/sqrt(numel(A0));
%     
%     c0 = arrayfun(@(t) t.endBuffer.sigma_r(mCh,k), tracks);
%     bufferBGMean(k+1) = mean(c0);
%     bufferBGStd(k+1) = std(c0)/sqrt(numel(c0));
% end
bufferMean = mean(lftData.endBuffer_Ia(:,:,mCh),1);
bufferStd = std(lftData.endBuffer_Ia(:,:,mCh),[],1) / sqrt(nIa);
% last frame of track
bufferBGMat = zeros(nIa,5);
for k = 1:nIa
    A0(k) = lftData.intMat_Ia(k,trackLengths_f(k),mCh);
    c0(k) = lftData.sigma_r_Ia(k,trackLengths_f(k)+nb,mCh);
    bufferBGMat(k,:) = lftData.sigma_r_Ia(k,trackLengths_f(k)+nb+(1:nb),mCh);    
end
bufferBGMean = mean(bufferBGMat,1);
bufferBGStd = std(bufferBGMat,[],1) / sqrt(nIa);

bufferMean = [mean(A0) bufferMean];
bufferStd = [std(A0)/sqrt(nIa) bufferStd];
bufferBGMean = [mean(c0) bufferBGMean];
bufferBGStd = [std(c0)/sqrt(nIa) bufferBGStd];
  
xlabels = [{'Last frame'} arrayfun(@(i) ['+' num2str(i*data.framerate) ' s'], 1:5, 'UniformOutput', false)];
figure;
barplot2([kLevel*bufferBGMean; bufferMean]', [kLevel*bufferBGStd; bufferStd]', 'ErrorbarPosition', 'both', ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1, 'YLim', YLim,...
    'XLabels', xlabels,...
    'YLabel', 'Mean intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'NorthOutside');
set(hl, 'Box', 'off');

if ip.Results.Print 
    print('-depsc2', 'EndBufferIntensities.eps');
end

%%

f1mean = zeros(nc,1);
f1std = zeros(nc,1);
f1bgmean = zeros(nc,1);
f1bgstd = zeros(nc,1);

intRatioPercentiles = zeros(5,nc);

bufferMean = zeros(nc,1);
bufferStd = zeros(nc,1);
bufferBGMean = zeros(nc,1);
bufferBGStd = zeros(nc,1);

for k = 1:nc
    cIdx = cohortBounds(k)<=trackLengths_f & trackLengths_f<cohortBounds(k+1);
    
    %A0 = arrayfun(@(t) t.A(mCh,1), tracks(cIdx));
    A0 = lftData.intMat_Ia(cIdx,1,mCh);
    f1mean(k) = mean(A0);
    f1std(k) = std(A0)/sqrt(nIa);
    %c0 = arrayfun(@(t) t.sigma_r(mCh,1), tracks(cIdx));
    c0 = lftData.sigma_r_Ia(cIdx,nb+1,mCh);
    f1bgmean(k) = mean(c0);
    f1bgstd(k) = std(c0)/sqrt(nIa);
    
    %maxInts = arrayfun(@(t) max(t.A(mCh,:)), tracks(cIdx));
    maxInts = nanmax(lftData.intMat_Ia(cIdx,:),[],2);
    % percentiles of the max. int. to 1st frame int. ratio
    intRatioPercentiles(:,k) = prctile(maxInts ./ A0, [50 25 75 2.5 97.5]);
    
    % mean start buffer int
    %A0 = arrayfun(@(t) t.startBuffer.A(mCh,end), tracks(cIdx));
    A0 = lftData.startBuffer_Ia(cIdx,end,mCh);
    bufferMean(k) = mean(A0);
    bufferStd(k) = std(A0)/sqrt(numel(A0));
    %c0 = arrayfun(@(t) t.startBuffer.sigma_r(mCh,end), tracks(cIdx));
    c0 = lftData.sigma_r_Ia(cIdx,nb+1,mCh);
    bufferBGMean(k) = mean(c0);
    bufferBGStd(k) = std(c0)/sqrt(numel(c0));
end



figure;
barplot2([kLevel*f1bgmean f1mean], [kLevel*f1bgstd f1std], 'ErrorbarPosition', 'both', ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1, 'YLim', YLim,...
    'XLabels', cohortLabels,...
    'XLabel', 'Lifetime cohort', 'YLabel', '1^{st} frame intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'NorthOutside');
set(hl, 'Box', 'off');
if ip.Results.Print 
    print('-depsc2', '1stFrameIntensitiesPerCohort.eps');
end


figure;
barplot2([kLevel*bufferBGMean bufferMean], [kLevel*bufferBGStd bufferStd], 'ErrorbarPosition', 'both', ...
    'FaceColor', [0.8 0.8 0.8; 0.7 0.9 1], 'EdgeColor', [0.6 0.6 0.6; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1, 'YLim', YLim,...
    'XLabels', cohortLabels,...
    'XLabel', 'Lifetime cohort', 'YLabel', 'Last start buffer frame intensity (A.U.)');
hl = legend('Background noise threshold', 'Intensity', 'Location', 'NorthOutside');
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

