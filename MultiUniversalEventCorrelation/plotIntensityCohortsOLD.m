function plotIntensityCohorts(data, intRes, condition)

if nargin<3 || isempty(condition)
    condition = '';
else
    condition = [', ' condition];
end

nMovies = length(data);
nCohorts = length(intRes);

colorV = ones(nCohorts,3);
colorV(:,1) = (nCohorts:-1:1)/nCohorts;
colorV = hsv2rgb(colorV);

figure;
for i = 1:nMovies;
    for c = 1:nCohorts
        % errorbars: standard error of the mean (SEM)
        errorbar(intRes(c).t(1,:), intRes(c).cIntensity(i,:), intRes(c).cIntensityStd(i,:)/sqrt(intRes(c).cohortSize(i)), 'Color', colorV(c,:));
        hold on;
    end    
end
set(gca, 'FontName', 'Helvetica', 'FontSize', 12, 'LineWidth', 1.5);
xlabel('time [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('intensity (above background) [A.U.]', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Cohorts for all movies' condition], 'FontName', 'Helvetica', 'FontSize', 14);


% plot mean w/ standard error of mean for all movies combined
figure;
for c = 1:nCohorts
    h = errorbar(intRes(c).framerate*intRes(c).tvec, intRes(c).intensityMean, intRes(c).intensitySEM, 'k-');
    h = get(h, 'Children');
    set(h(1), 'LineWidth', 2, 'Color', colorV(c,:));    
    hold on;
end
set(gca, 'FontName', 'Helvetica', 'FontSize', 12, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Intensity (A.U.)', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Averaged cohorts' condition], 'FontName', 'Helvetica', 'FontSize', 14);
