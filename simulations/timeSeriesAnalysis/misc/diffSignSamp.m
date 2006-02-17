
numTrials = 20000;
trajLength = [100:50:400]';
numTraj = 10;

numSamps = length(trajLength);

meanNumGreater = zeros(numSamps,1);
varNumGreater = zeros(numSamps,1);

for samp = 1:numSamps

    numGreater = zeros(numTrials,1);

    for trial = 1:numTrials

        for i=1:numTraj

            x1 = randn(trajLength(samp)+1,1);
            xDiff = x1(2:end) - x1(1:end-1);
            numGreater(trial) = numGreater(trial) + length(find(xDiff>0));

        end

    end

    meanNumGreater(samp) = mean(numGreater);
    varNumGreater(samp) = std(numGreater)^2;

end