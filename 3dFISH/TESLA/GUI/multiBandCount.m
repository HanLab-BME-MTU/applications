function bandStat = multiBandCount(bandStat)
%MULTIBANDCOUNT Selects the bands with high intensity and assign them
%multiple counting numbers
%   Detailed explanation goes here

yPos = zeros(numel(bandStat), 1);
for i = 1:numel(bandStat)
    bandStat(i).yPos = bandStat(i).bandPos(2);
    yPos(i) = bandStat(i).yPos;
end

interval = round((max(yPos) - min(yPos))/20);
for bandIndex = 1:numel(bandStat)
    botLim = min(yPos(bandIndex) + interval, max(yPos));
    topLim = max(yPos(bandIndex) - interval, min(yPos));
    subBandStat = bandStat([bandStat(:).yPos] < botLim & [bandStat(:).yPos] > topLim);
    % Choice between median intensity and average intensity
    medIntensity = median([subBandStat(:).bandIntensity]);
    bandStat(bandIndex).count = max(1, round(bandStat(bandIndex).bandIntensity/medIntensity));
end

end