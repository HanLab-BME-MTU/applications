function [bandStat, shortThreshPos, marker] = bandStatCalcPWL(handles, fIM, shortThresh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

bandMap = handles.bandMap;
ladder = handles.ladder;
halfBandRange = handles.halfBandRange;
marker.num = handles.markerNum;
marker.pos = handles.markerPos;
marker.size = handles.markerSize;

bandStat = struct('index', [], 'bandPos', [], 'bandSize', [], 'bandIntensity', [], 'count', []);

% Linear regression y = ax+b for any two adjacent marker pair (x:marker.pos, y:marker.size)
G(:,1) = marker.pos;
G(:,2) = ones(marker.num, 1);
para = struct('slope', [], 'intercept', []);
for num = 1:marker.num-1
    regPara = G(num:num+1,:)\marker.size(num:num+1);
    para(num).slope = regPara(1);
    para(num).intercept = regPara(2);
end

% Threshold position annotation
if shortThresh >= marker.size(1)
    shortThreshPos = (shortThresh-para(1).intercept)/para(1).slope;
else if shortThresh < marker.size(marker.num)
        shortThreshPos = (shortThresh-para(marker.num).intercept)/para(marker.num).slope;
    else
        for num = 1:marker.num-1
            if shortThresh < marker.size(num) && shortThresh >= marker.size(num+1)
                shortThreshPos = (shortThresh-para(num).intercept)/para(num).slope;
            end
        end
    end
end

% Band size annotation
bandIndex = 0;
for q=1:size(fIM,2)
    for p=1:size(fIM,1)
        if bandMap(p,q) == 1 && (q < ladder.left || q > ladder.right)
            bandIndex = bandIndex + 1;
            bandStat(bandIndex).index = bandIndex;
            % Bands positions are recorded in cartesian coordinate
            bandStat(bandIndex).bandPos = [q, p];
            bandStat(bandIndex).count = 1;
            
            if p <= marker.pos(1)
                bandStat(bandIndex).bandSize = para(1).slope*p+para(1).intercept;
            else if p > marker.pos(marker.num)
                    bandStat(bandIndex).bandSize = para(marker.num-1).slope*p+para(marker.num-1).intercept;
                else
                    for num = 1:marker.num-1
                        if p > marker.pos(num) && p <= marker.pos(num+1)
                            bandStat(bandIndex).bandSize = para(num).slope*p+para(num).intercept;
                        end
                    end
                end
            end
            
            % Calculate average intensity surrounding the seed
            halfBandHeight = round(0.2*halfBandRange);
            bandIntensity = 0;
            bandLeftBoundary = max(0, q - halfBandRange);
            bandRightBoundary = min(size(fIM,2), q + halfBandRange);
            bandTopBoundary = max(0, p - halfBandHeight);
            bandBotBoundary = min (size(fIM,1), p + halfBandHeight);
            count = 0;
            for horiRange = bandLeftBoundary:bandRightBoundary
                for vertiRange = bandTopBoundary:bandBotBoundary
                    bandIntensity = bandIntensity + fIM(vertiRange, horiRange);
                    count = count + 1;
                end
            end
            
            bandStat(bandIndex).bandIntensity = bandIntensity/count;
        end
    end
end

for trivialBand = 1:size(bandStat,2)
    if bandStat(trivialBand).bandSize <= 0
        bandStat(trivialBand).bandSize = 0.00001;
    end
end

end