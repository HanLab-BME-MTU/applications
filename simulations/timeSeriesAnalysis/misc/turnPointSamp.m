
numTrials = 20000;
trajLength = [10000:1000:15000]';
numTraj = 1;

numSamps = length(trajLength);

meanNumTurn = zeros(numSamps,1);
varNumTurn = zeros(numSamps,1);

for samp = 1:numSamps

    numTurn = zeros(numTrials,1);

    for trial = 1:numTrials

        for i=1:numTraj

            x1 = randn(trajLength(samp),1);
            xDiff = (x1(2:end-1)-x1(1:end-2)).*(x1(2:end-1)-x1(3:end));
            numTurn(trial) = numTurn(trial) + length(find(xDiff>0));

        end

    end

    meanNumTurn(samp) = mean(numTurn);
    varNumTurn(samp) = std(numTurn)^2;

end

% x2 = randn(trajLength(samp),1);
% x3 = randn(trajLength(samp),1);
% xDiff = (x1-x2).*(x1-x3);
