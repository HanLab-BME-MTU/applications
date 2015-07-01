function corr = trajCorrelation(samp1,samp2,time,avgWin,lagT)

%% Setup
samp1 = squeeze(samp1);
samp2 = squeeze(samp2);

if size(samp1) == size(samp2')
    samp2 = samp2';
end

[numFrames, dim] = size(samp1);
if numFrames < dim
    samp1 = samp1';
    samp2 = samp2';
    [numFrames, ~] = size(samp1);
end

totalTime = time(end) - time(1) + 1;

%% Calculate Correlation Coefficient

if (totalTime == 1)
    corr = 0;
    startT = time - round((avgWin+lagT)/2);
    if (startT + avgWin + lagT) <= numFrames
        samp1Vec = zeros(avgWin,2);
        samp2Vec = zeros(avgWin,2);
        
        for i = 1:avgWin
            samp1Vec(i,:) = samp1(i+lagT+startT,:)-samp1(i+startT,:);
            samp2Vec(i,:) = samp2(i+lagT+startT,:)-samp2(i+startT,:);
        end
        dotProducts = sum(samp1Vec.*samp2Vec,2);
        corr = sum(dotProducts)/sqrt(sum(sum(samp1Vec.^2))*sum(sum(samp2Vec.^2)));
    end
else
    corr = zeros(length(time),1);
    if ~isempty(time)
        samp1Vec = zeros(avgWin+totalTime-1,2);
        samp2Vec = zeros(avgWin+totalTime-1,2);
        startT = time - round((avgWin+lagT)/2);
       
        for i = 1:(avgWin+totalTime-1)
            samp1Vec(i,:) = samp1(i+lagT+startT(1),:)-samp1(i+startT(1),:);
            samp2Vec(i,:) = samp2(i+lagT+startT(1),:)-samp2(i+startT(1),:);
        end
        dotProducts = sum(samp1Vec.*samp2Vec,2);
        
        for t = 1:length(startT)
            if (startT(t) + avgWin + lagT) <= numFrames
                corr(t) = sum(dotProducts(t:(t+avgWin-1)))/sqrt(sum(sum(samp1Vec(t:(t+avgWin-1),:).^2))*sum(sum(samp2Vec(t:(t+avgWin-1),:).^2)));
            end
        end
    end
end